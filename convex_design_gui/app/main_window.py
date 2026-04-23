import json
import math
import os
import random
import re
import shutil
import signal
import subprocess
import sys
import time
import csv
from pathlib import Path
from typing import Dict, Iterable, List, Optional, Sequence, Tuple

if os.name == "posix":
    import select

from PySide6.QtCore import QObject, Qt, QThread, Signal
from PySide6.QtWidgets import (
    QApplication,
    QComboBox,
    QDialog,
    QDialogButtonBox,
    QFileDialog,
    QFormLayout,
    QGridLayout,
    QGroupBox,
    QHBoxLayout,
    QHeaderView,
    QLabel,
    QLineEdit,
    QMainWindow,
    QMessageBox,
    QPlainTextEdit,
    QProgressBar,
    QPushButton,
    QSpinBox,
    QTableWidget,
    QTableWidgetItem,
    QVBoxLayout,
    QWidget,
)

from Bio.PDB import MMCIFParser, PDBIO, PDBParser, Select, is_aa


PROJECT_ROOT = Path(__file__).resolve().parents[1]
SETTINGS_PATH = PROJECT_ROOT / "resources" / "settings.json"
SCRIPT_COLUMNS = [
    "input file",
    "target chain",
    "sphere convexity",
    "sphere radius",
    "surface shape",
    "library length",
    "hotspot residues",
    "number of designs",
    "gpu",
]
SCRIPT_FIELD_LABELS = [
    ("input file", "File"),
    ("target chain", "Chain"),
    ("sphere convexity", "Convexity"),
    ("sphere radius", "Radius"),
    ("surface shape", "Shape"),
    ("library length", "Lib Length"),
    ("hotspot residues", "Hotspots"),
    ("number of designs", "Designs"),
    ("gpu", "GPU"),
]
OUTPUT_COLUMNS = [
    "pdb file",
    "sphere convexity",
    "convexity difference",
    "shape match",
    "length",
    "status",
]
OUTPUT_HEADER_LABELS = {
    "pdb file": "PDB",
    "sphere convexity": "Convexity",
    "convexity difference": "Diff",
    "shape match": "Match",
    "length": "Len",
    "status": "Status",
}
AA3_TO_1 = {
    "ALA": "A",
    "ARG": "R",
    "ASN": "N",
    "ASP": "D",
    "CYS": "C",
    "GLN": "Q",
    "GLU": "E",
    "GLY": "G",
    "HIS": "H",
    "ILE": "I",
    "LEU": "L",
    "LYS": "K",
    "MET": "M",
    "PHE": "F",
    "PRO": "P",
    "SER": "S",
    "THR": "T",
    "TRP": "W",
    "TYR": "Y",
    "VAL": "V",
}


class RunSignals(QObject):
    log = Signal(str)
    error = Signal(str)
    status = Signal(str)
    finished = Signal(int)
    output_count = Signal(int)
    output_pdb = Signal(dict)
    progress = Signal(int)


class ScriptRunWorker(QObject):
    def __init__(
        self,
        script_path: Path,
        output_dir: Path,
        expected_designs: int,
        target_chain: str,
        target_sequence: str,
        target_length: int,
    ) -> None:
        super().__init__()
        self.script_path = script_path
        self.output_dir = output_dir
        self.expected_designs = max(1, int(expected_designs))
        self.target_chain = target_chain
        self.target_sequence = target_sequence
        self.target_length = target_length
        self.signals = RunSignals()
        self._process = None  # type: Optional[subprocess.Popen]
        self._cancel_requested = False
        self._seen_outputs = set()  # type: set
        self._output_size_cache = {}  # type: Dict[str, int]
        self._expected_count_seen_at = None  # type: Optional[float]

    def request_cancel(self) -> None:
        self._cancel_requested = True
        if self._process is not None and self._process.poll() is None:
            self.signals.log.emit("[INFO] Terminating active job.")
            try:
                self._process.terminate()
            except Exception:
                pass

    def run(self) -> None:
        self.signals.status.emit("running")
        try:
            command = ["bash", str(self.script_path)]
            self.signals.log.emit("[INFO] Running command: " + " ".join(command))
            self._process = subprocess.Popen(
                command,
                cwd=str(self.script_path.parent),
                stdout=subprocess.PIPE,
                stderr=subprocess.STDOUT,
                text=True,
                bufsize=1,
                start_new_session=True,
            )
            completed_by_output = self._monitor_process_until_done()
            return_code = self._wait_for_process()
            if completed_by_output and return_code != 0:
                return_code = 0
            self._scan_outputs(status="Completed", force=True)
            if self._cancel_requested:
                self.signals.status.emit("cancelled")
            elif return_code == 0:
                self.signals.status.emit("Completed")
            else:
                self.signals.status.emit("failed")
            self.signals.output_count.emit(self._count_output_pdbs())
            self.signals.finished.emit(int(return_code))
        except Exception as exc:
            self.signals.error.emit(str(exc))
            self.signals.status.emit("failed")
            self.signals.output_count.emit(self._count_output_pdbs())
            self.signals.finished.emit(1)

    def _monitor_process_until_done(self) -> bool:
        completed_by_output = False
        output_count = 0
        stable_count = 0
        last_scan_at = 0.0
        scan_interval = 1.0
        while self._process is not None:
            read_output = self._read_available_output()
            now = time.monotonic()
            process_done = self._process.poll() is not None
            if process_done or now - last_scan_at >= scan_interval:
                output_count, stable_count = self._scan_outputs(status="Completed", force=False)
                self.signals.output_count.emit(output_count)
                last_scan_at = now

            if self._cancel_requested:
                self._terminate_process()
                break

            if output_count >= self.expected_designs:
                if self._expected_count_seen_at is None:
                    self._expected_count_seen_at = now
                elif now - self._expected_count_seen_at >= 2.0 and stable_count >= self.expected_designs:
                    completed_by_output = True
                    self.signals.log.emit(
                        f"[INFO] Expected output count reached ({output_count}/{self.expected_designs}). "
                        "Stopping RFdiffusion shutdown wait and finalizing results."
                    )
                    self._scan_outputs(status="Completed", force=True)
                    self._terminate_process()
                    break
            else:
                self._expected_count_seen_at = None

            if self._process.poll() is not None:
                self._drain_remaining_output()
                break

            if not read_output:
                time.sleep(0.05)
        return completed_by_output

    def _read_available_output(self) -> bool:
        if self._process is None or self._process.stdout is None:
            return False
        stdout = self._process.stdout
        read_any = False
        if os.name == "posix":
            while True:
                ready, _, _ = select.select([stdout], [], [], 0)
                if not ready:
                    break
                raw_line = stdout.readline()
                if not raw_line:
                    break
                line = raw_line.rstrip()
                if line:
                    self.signals.log.emit(line)
                    read_any = True
            return read_any
        elif self._process.poll() is None:
            return False
        raw_line = stdout.readline()
        if raw_line:
            line = raw_line.rstrip()
            if line:
                self.signals.log.emit(line)
                read_any = True
        return read_any

    def _drain_remaining_output(self) -> None:
        if self._process is None or self._process.stdout is None:
            return
        for raw_line in self._process.stdout:
            line = raw_line.rstrip()
            if line:
                self.signals.log.emit(line)

    def _wait_for_process(self) -> int:
        if self._process is None:
            return 1
        try:
            return int(self._process.wait(timeout=10))
        except subprocess.TimeoutExpired:
            self._terminate_process(force=True)
            return int(self._process.wait(timeout=10))

    def _terminate_process(self, force: bool = False) -> None:
        if self._process is None or self._process.poll() is not None:
            return
        try:
            if os.name == "posix":
                sig = signal.SIGKILL if force else signal.SIGTERM
                os.killpg(os.getpgid(self._process.pid), sig)
            elif force:
                self._process.kill()
            else:
                self._process.terminate()
        except Exception:
            try:
                self._process.kill() if force else self._process.terminate()
            except Exception:
                pass

    def _count_output_pdbs(self) -> int:
        if not self.output_dir.exists():
            return 0
        return sum(1 for path in self._iter_output_pdbs())

    def _iter_output_pdbs(self):
        if not self.output_dir.exists():
            return []
        return [
            path
            for path in self.output_dir.rglob("*.pdb")
            if path.is_file() and "traj" not in path.relative_to(self.output_dir).parts
        ]

    def _scan_outputs(self, status: str, force: bool = False) -> Tuple[int, int]:
        output_paths = sorted(self._iter_output_pdbs())
        stable_count = 0
        for path in output_paths:
            key = str(path)

            try:
                current_size = path.stat().st_size
            except OSError:
                continue
            previous_size = self._output_size_cache.get(key)
            self._output_size_cache[key] = current_size
            is_stable = previous_size == current_size
            if is_stable:
                stable_count += 1

            if key in self._seen_outputs:
                continue

            if not force and previous_size != current_size:
                continue

            metrics = calculate_output_pdb_metrics(
                path,
                self.target_chain,
                target_sequence=self.target_sequence,
                target_length=self.target_length,
            )
            metrics_ready = bool(metrics.get("sphere convexity") not in ("", None) and metrics.get("length") not in ("", None))
            if not metrics_ready and not force:
                continue

            self._seen_outputs.add(key)
            self.signals.output_pdb.emit(
                {
                    "pdb file": path.name,
                    "sphere convexity": metrics.get("sphere convexity", ""),
                    "length": metrics.get("length", ""),
                    "target chain": metrics.get("target chain", ""),
                    "binder chain": metrics.get("binder chain", ""),
                    "path": str(path),
                    "status": status,
                }
            )
        progress = int(min(100, round((len(output_paths) / self.expected_designs) * 100)))
        self.signals.progress.emit(progress)
        return len(output_paths), stable_count


def load_settings() -> Dict[str, str]:
    defaults = {
        "hcs_dir": "/data/5HCS",
        "hcs_extracted_dir": "/data/5HCS/5HCS_extracted",
    }
    if not SETTINGS_PATH.exists():
        return defaults
    try:
        loaded = json.loads(SETTINGS_PATH.read_text(encoding="utf-8"))
    except Exception:
        return defaults
    if not isinstance(loaded, dict):
        return defaults
    for key in defaults:
        value = loaded.get(key, defaults[key])
        defaults[key] = str(value).strip() or defaults[key]
    return defaults


def save_settings(settings: Dict[str, str]) -> None:
    SETTINGS_PATH.parent.mkdir(parents=True, exist_ok=True)
    SETTINGS_PATH.write_text(json.dumps(settings, indent=2), encoding="utf-8")


def calculate_paper_style_convexity(
    input_pdb: Path,
    chain_id: str,
    hotspot_residues: Sequence[str],
) -> Dict[str, object]:
    structure = _load_structure(input_pdb)
    model = structure[0]
    if chain_id not in model:
        available = ", ".join(chain.id for chain in model)
        raise ValueError(f"Chain '{chain_id}' was not found in {input_pdb.name}. Available chains: {available or '(none)'}")

    hotspot_numbers = _hotspot_numbers_for_chain(hotspot_residues, chain_id)
    chain = model[chain_id]
    chain_heavy_coords = _collect_heavy_atom_coords(chain, residue_numbers=set())
    patch_heavy_coords = _collect_heavy_atom_coords(chain, residue_numbers=hotspot_numbers)
    if len(patch_heavy_coords) < 4:
        raise ValueError("At least four hotspot heavy atoms are required for paper-style convexity.")

    return _calculate_signed_convexity_from_interface(
        protein_heavy_coords=chain_heavy_coords,
        interface_heavy_coords=patch_heavy_coords,
    )


def _calculate_signed_convexity_from_interface(
    protein_heavy_coords: Sequence[Tuple[float, float, float]],
    interface_heavy_coords: Sequence[Tuple[float, float, float]],
) -> Dict[str, object]:
    if len(interface_heavy_coords) < 4 or len(protein_heavy_coords) < 4:
        return {
            "sphere radius": "",
            "sphere convexity": "",
            "surface shape": "",
        }

    sphere = _ransac_sphere_fit(
        coords=interface_heavy_coords,
        threshold=1.0,
        max_iterations=100000,
        seed=17,
    )
    if sphere is None:
        return {
            "sphere radius": "",
            "sphere convexity": "",
            "surface shape": "",
        }

    sphere_center, radius = sphere
    protein_center = _centroid(protein_heavy_coords)
    interface_center = _centroid(interface_heavy_coords)
    protein_vector = _vector(interface_center, protein_center)
    fitted_vector = _vector(interface_center, sphere_center)
    direction_score = _dot(protein_vector, fitted_vector)
    magnitude = 1.0 / radius if radius > 0.0 else ""
    signed_convexity = -magnitude if direction_score >= 0.0 else magnitude
    return {
        "sphere radius": radius,
        "sphere convexity": signed_convexity,
        "surface shape": "convex" if direction_score >= 0.0 else "concave",
    }


def calculate_output_pdb_metrics(
    output_pdb: Path,
    target_chain: str,
    target_sequence: str = "",
    target_length: int = 0,
) -> Dict[str, object]:
    try:
        structure = _load_structure(output_pdb)
        model = structure[0]
        target = _select_output_target_chain(
            model=model,
            preferred_chain_id=target_chain,
            target_sequence=target_sequence,
            target_length=target_length,
        )
        if target is None:
            return {"sphere convexity": "", "length": "", "target chain": "", "binder chain": ""}

        target_coords = _collect_heavy_atom_coords(target, residue_numbers=set())
        candidate_chains = [chain for chain in model if chain.id != target.id]
        if not candidate_chains or not target_coords:
            return {"sphere convexity": "", "length": "", "target chain": target.id, "binder chain": ""}

        chain = max(candidate_chains, key=lambda candidate: len(_collect_heavy_atom_coords(candidate, residue_numbers=set())))
        binder_coords = _collect_heavy_atom_coords(chain, residue_numbers=set())
        interface_coords = _collect_interface_heavy_atom_coords(
            query_chain=chain,
            partner_coords=target_coords,
            cutoff=5.0,
        )
        residue_numbers = set()
        for residue in chain:
            if residue.id[0] == " " and is_aa(residue, standard=False) and "CA" in residue:
                residue_numbers.add(int(residue.id[1]))
        metrics = _calculate_signed_convexity_from_interface(
            protein_heavy_coords=binder_coords,
            interface_heavy_coords=interface_coords,
        )
        return {
            "sphere convexity": metrics.get("sphere convexity", ""),
            "length": len(residue_numbers),
            "target chain": target.id,
            "binder chain": chain.id,
        }
    except Exception:
        return {"sphere convexity": "", "length": "", "target chain": "", "binder chain": ""}


def _select_output_target_chain(
    model,
    preferred_chain_id: str,
    target_sequence: str,
    target_length: int,
):
    chains = [chain for chain in model if _count_chain_residues(chain) > 0]
    if not chains:
        return None

    if target_sequence:
        return max(
            chains,
            key=lambda chain: (
                _sequence_similarity(_chain_sequence(chain), target_sequence),
                -abs(_count_chain_residues(chain) - len(target_sequence)),
            ),
        )

    if preferred_chain_id in model:
        return model[preferred_chain_id]

    if target_length > 0:
        return min(chains, key=lambda chain: abs(_count_chain_residues(chain) - target_length))
    return chains[0]


def _sequence_similarity(first: str, second: str) -> float:
    if not first or not second:
        return 0.0
    length = min(len(first), len(second))
    matches = sum(1 for index in range(length) if first[index] == second[index])
    length_penalty = abs(len(first) - len(second)) / max(len(first), len(second), 1)
    return (matches / length) - length_penalty


def _load_structure(input_pdb: Path):
    suffix = input_pdb.suffix.lower()
    if suffix == ".pdb":
        parser = PDBParser(QUIET=True)
    elif suffix in {".cif", ".mmcif"}:
        parser = MMCIFParser(QUIET=True)
    else:
        raise ValueError(f"Unsupported structure format: {input_pdb.suffix}")
    return parser.get_structure(input_pdb.stem, str(input_pdb))


def _chain_signature(input_pdb: Path, chain_id: str) -> Tuple[int, str]:
    structure = _load_structure(input_pdb)
    model = structure[0]
    if chain_id not in model:
        available = ", ".join(chain.id for chain in model)
        raise ValueError(f"Chain '{chain_id}' was not found in {input_pdb.name}. Available chains: {available or '(none)'}")
    chain = model[chain_id]
    return _count_chain_residues(chain), _chain_sequence(chain)


def _count_chain_residues(chain) -> int:
    return sum(1 for residue in chain if residue.id[0] == " " and is_aa(residue, standard=False) and "CA" in residue)


def _chain_sequence(chain) -> str:
    letters = []
    for residue in chain:
        if residue.id[0] != " " or not is_aa(residue, standard=False) or "CA" not in residue:
            continue
        letters.append(AA3_TO_1.get(residue.get_resname().upper(), "X"))
    return "".join(letters)


class _SingleChainSelect(Select):
    def __init__(self, chain_id: str) -> None:
        self.chain_id = chain_id

    def accept_model(self, model) -> int:
        return 1 if model.id == 0 else 0

    def accept_chain(self, chain) -> int:
        return 1 if chain.id == self.chain_id else 0


def _write_target_chain_pdb(input_pdb: Path, chain_id: str, output_pdb: Path) -> None:
    structure = _load_structure(input_pdb)
    model = structure[0]
    if chain_id not in model:
        available = ", ".join(chain.id for chain in model)
        raise ValueError(f"Chain '{chain_id}' was not found in {input_pdb.name}. Available chains: {available or '(none)'}")

    output_pdb.parent.mkdir(parents=True, exist_ok=True)
    io = PDBIO()
    io.set_structure(structure)
    io.save(str(output_pdb), select=_SingleChainSelect(chain_id))

    if not output_pdb.is_file() or output_pdb.stat().st_size == 0:
        raise ValueError(f"Failed to write target chain PDB: {output_pdb}")


def _hotspot_numbers_for_chain(hotspot_residues: Sequence[str], chain_id: str) -> set:
    numbers = set()
    for hotspot in hotspot_residues:
        match = re.fullmatch(r"([A-Za-z])(\d+)", hotspot)
        if match and match.group(1).upper() == chain_id:
            numbers.add(int(match.group(2)))
    return numbers


def _collect_heavy_atom_coords(chain, residue_numbers: set) -> List[Tuple[float, float, float]]:
    coords = []
    for residue in chain:
        if residue.id[0] != " " or not is_aa(residue, standard=False):
            continue
        residue_number = int(residue.id[1])
        if residue_numbers and residue_number not in residue_numbers:
            continue
        for atom in residue:
            element = str(getattr(atom, "element", "")).strip().upper()
            if element == "H":
                continue
            coord = atom.get_coord()
            coords.append((float(coord[0]), float(coord[1]), float(coord[2])))
    return coords


def _collect_interface_heavy_atom_coords(
    query_chain,
    partner_coords: Sequence[Tuple[float, float, float]],
    cutoff: float,
) -> List[Tuple[float, float, float]]:
    cutoff_squared = cutoff * cutoff
    interface_coords = []
    for residue in query_chain:
        if residue.id[0] != " " or not is_aa(residue, standard=False):
            continue
        for atom in residue:
            element = str(getattr(atom, "element", "")).strip().upper()
            if element == "H":
                continue
            coord_array = atom.get_coord()
            coord = (float(coord_array[0]), float(coord_array[1]), float(coord_array[2]))
            for partner_coord in partner_coords:
                squared_distance = (
                    (coord[0] - partner_coord[0]) ** 2
                    + (coord[1] - partner_coord[1]) ** 2
                    + (coord[2] - partner_coord[2]) ** 2
                )
                if squared_distance <= cutoff_squared:
                    interface_coords.append(coord)
                    break
    return interface_coords


def _ransac_sphere_fit(
    coords: Sequence[Tuple[float, float, float]],
    threshold: float,
    max_iterations: int,
    seed: int,
) -> Optional[Tuple[Tuple[float, float, float], float]]:
    if len(coords) < 4:
        return None
    rng = random.Random(seed)
    best_center = None
    best_radius = 0.0
    best_inliers = -1
    coord_list = list(coords)
    iterations = min(max_iterations, max(200, len(coord_list) * 20))

    for _ in range(iterations):
        sphere = _sphere_from_four_points(rng.sample(coord_list, 4))
        if sphere is None:
            continue
        center, radius = sphere
        if radius <= 0.0 or not math.isfinite(radius):
            continue
        inliers = sum(1 for coord in coord_list if abs(_distance(coord, center) - radius) <= threshold)
        if inliers > best_inliers:
            best_center = center
            best_radius = radius
            best_inliers = inliers

    if best_center is None:
        return None
    return best_center, best_radius


def _sphere_from_four_points(
    points: Sequence[Tuple[float, float, float]]
) -> Optional[Tuple[Tuple[float, float, float], float]]:
    matrix = []
    values = []
    for x, y, z in points:
        matrix.append([x, y, z, 1.0])
        values.append(-(x * x + y * y + z * z))
    solution = _solve_linear_system(matrix, values)
    if solution is None:
        return None
    a_value, b_value, c_value, d_value = solution
    center = (-a_value / 2.0, -b_value / 2.0, -c_value / 2.0)
    radius_squared = center[0] ** 2 + center[1] ** 2 + center[2] ** 2 - d_value
    if radius_squared <= 0.0:
        return None
    return center, math.sqrt(radius_squared)


def _solve_linear_system(matrix: List[List[float]], values: List[float]) -> Optional[List[float]]:
    size = len(values)
    augmented = [row[:] + [values[index]] for index, row in enumerate(matrix)]
    for col in range(size):
        pivot = max(range(col, size), key=lambda row_index: abs(augmented[row_index][col]))
        if abs(augmented[pivot][col]) < 1e-9:
            return None
        augmented[col], augmented[pivot] = augmented[pivot], augmented[col]
        pivot_value = augmented[col][col]
        for item_index in range(col, size + 1):
            augmented[col][item_index] /= pivot_value
        for row_index in range(size):
            if row_index == col:
                continue
            factor = augmented[row_index][col]
            for item_index in range(col, size + 1):
                augmented[row_index][item_index] -= factor * augmented[col][item_index]
    return [augmented[row_index][size] for row_index in range(size)]


def _centroid(coords: Iterable[Tuple[float, float, float]]) -> Tuple[float, float, float]:
    values = list(coords)
    if not values:
        return (0.0, 0.0, 0.0)
    return (
        sum(coord[0] for coord in values) / len(values),
        sum(coord[1] for coord in values) / len(values),
        sum(coord[2] for coord in values) / len(values),
    )


def _distance(first: Tuple[float, float, float], second: Tuple[float, float, float]) -> float:
    return math.sqrt(
        (first[0] - second[0]) ** 2
        + (first[1] - second[1]) ** 2
        + (first[2] - second[2]) ** 2
    )


def _vector(first: Tuple[float, float, float], second: Tuple[float, float, float]) -> Tuple[float, float, float]:
    return (second[0] - first[0], second[1] - first[1], second[2] - first[2])


def _dot(first: Tuple[float, float, float], second: Tuple[float, float, float]) -> float:
    return first[0] * second[0] + first[1] * second[1] + first[2] * second[2]


class SettingsDialog(QDialog):
    def __init__(self, settings: Dict[str, str], parent: Optional[QWidget] = None) -> None:
        super().__init__(parent)
        self.setWindowTitle("Settings")
        self.resize(760, 180)
        self._build_ui(settings)

    def _build_ui(self, settings: Dict[str, str]) -> None:
        layout = QVBoxLayout(self)
        form = QFormLayout()
        form.setContentsMargins(8, 8, 8, 8)
        form.setSpacing(10)

        self.hcs_dir_edit = self._create_directory_row(
            form=form,
            label="HCS_DIR",
            value=settings.get("hcs_dir", ""),
            title="Select HCS Directory",
        )
        self.hcs_extracted_dir_edit = self._create_directory_row(
            form=form,
            label="HCS_EXTRACTED_DIR",
            value=settings.get("hcs_extracted_dir", ""),
            title="Select Extracted HCS Directory",
        )

        layout.addLayout(form)
        buttons = QDialogButtonBox(QDialogButtonBox.Save | QDialogButtonBox.Cancel)
        buttons.accepted.connect(self.accept)
        buttons.rejected.connect(self.reject)
        layout.addWidget(buttons)

    def _create_directory_row(self, form: QFormLayout, label: str, value: str, title: str) -> QLineEdit:
        row = QWidget()
        row_layout = QHBoxLayout(row)
        row_layout.setContentsMargins(0, 0, 0, 0)
        row_layout.setSpacing(6)

        line_edit = QLineEdit()
        line_edit.setText(value)
        line_edit.setPlaceholderText(title)

        browse_btn = QPushButton("Browse")
        browse_btn.setFixedWidth(90)

        def browse() -> None:
            selected_dir = QFileDialog.getExistingDirectory(self, title, line_edit.text().strip())
            if selected_dir:
                line_edit.setText(selected_dir)

        browse_btn.clicked.connect(browse)
        row_layout.addWidget(line_edit, 1)
        row_layout.addWidget(browse_btn)
        form.addRow(label, row)
        return line_edit

    def get_settings(self) -> Dict[str, str]:
        return {
            "hcs_dir": self.hcs_dir_edit.text().strip(),
            "hcs_extracted_dir": self.hcs_extracted_dir_edit.text().strip(),
        }


class MainWindow(QMainWindow):
    def __init__(self) -> None:
        super().__init__()
        self.settings = load_settings()
        self.run_thread = None  # type: Optional[QThread]
        self.run_worker = None  # type: Optional[ScriptRunWorker]
        self.current_input_convexity = 0.0
        self.current_top_n = 1
        self.current_selection_dir = None  # type: Optional[Path]
        self.output_rows = []  # type: List[Dict[str, object]]
        self.setWindowTitle("Convex Design GUI")
        self.resize(1240, 820)
        self._build_ui()
        self._fit_result_tables()

    def _build_ui(self) -> None:
        central = QWidget()
        root = QVBoxLayout(central)
        root.setContentsMargins(16, 16, 16, 16)
        root.setSpacing(10)

        header = QHBoxLayout()
        title = QLabel("Convex Design GUI")
        title.setStyleSheet("font-size: 22px; font-weight: 700;")
        header.addWidget(title)
        header.addStretch(1)

        self.settings_btn = QPushButton("Settings")
        self.settings_btn.setMinimumHeight(34)
        self.settings_btn.clicked.connect(self._open_settings_dialog)
        header.addWidget(self.settings_btn)
        root.addLayout(header)

        top_row = QHBoxLayout()
        top_row.setSpacing(10)
        top_row.addWidget(self._build_input_group(), 3)
        top_row.addWidget(self._build_run_group(), 2)
        root.addLayout(top_row)

        root.addWidget(self._build_input_result_group())
        root.addWidget(self._build_output_result_group())

        self.log_edit = QPlainTextEdit()
        self.log_edit.setReadOnly(True)
        self.log_edit.setMaximumHeight(230)
        root.addWidget(self.log_edit)

        self.progress_bar = QProgressBar()
        self.progress_bar.setRange(0, 100)
        self.progress_bar.setValue(0)
        root.addWidget(self.progress_bar)

        self.setCentralWidget(central)

    def _build_input_group(self) -> QGroupBox:
        group = QGroupBox("RFdiffusion Input")
        layout = QGridLayout(group)
        layout.setContentsMargins(10, 12, 10, 10)
        layout.setHorizontalSpacing(8)
        layout.setVerticalSpacing(8)

        layout.addWidget(QLabel("Input PDB"), 0, 0)
        self.input_pdb_edit = QLineEdit()
        self.input_pdb_edit.setPlaceholderText("Target-only PDB file")
        layout.addWidget(self.input_pdb_edit, 0, 1)

        input_browse_btn = QPushButton("Browse")
        input_browse_btn.clicked.connect(self._browse_input_pdb)
        layout.addWidget(input_browse_btn, 0, 2)

        layout.addWidget(QLabel("Target chain"), 1, 0)
        chain_row = QWidget()
        chain_row_layout = QHBoxLayout(chain_row)
        chain_row_layout.setContentsMargins(0, 0, 0, 0)
        chain_row_layout.setSpacing(4)

        self.target_chain_combo = QComboBox()
        self.target_chain_combo.addItems(["A", "B", "C", "D", "E", "Other"])
        self.target_chain_combo.setFixedWidth(86)
        self.target_chain_combo.currentTextChanged.connect(self._on_target_chain_changed)
        chain_row_layout.addWidget(self.target_chain_combo)

        self.target_chain_other_edit = QLineEdit()
        self.target_chain_other_edit.setPlaceholderText("Enter chain ID")
        self.target_chain_other_edit.setFixedWidth(110)
        self.target_chain_other_edit.setVisible(False)
        chain_row_layout.addWidget(self.target_chain_other_edit)
        chain_row_layout.addStretch(1)
        layout.addWidget(chain_row, 1, 1, 1, 2)

        layout.addWidget(QLabel("Library search length"), 2, 0)
        self.library_length_edit = QLineEdit("70-140")
        self.library_length_edit.setPlaceholderText("Example: 70-140")
        layout.addWidget(self.library_length_edit, 2, 1, 1, 2)

        layout.addWidget(QLabel("Hotspot residues"), 3, 0)
        self.hotspot_residue_edit = QLineEdit()
        self.hotspot_residue_edit.setPlaceholderText("Example: A96,A98,A106 or 96,98,106")
        layout.addWidget(self.hotspot_residue_edit, 3, 1, 1, 2)

        layout.addWidget(QLabel("Number of designs"), 4, 0)
        design_row = QWidget()
        design_row_layout = QHBoxLayout(design_row)
        design_row_layout.setContentsMargins(0, 0, 0, 0)
        design_row_layout.setSpacing(8)

        self.num_designs_spin = QSpinBox()
        self.num_designs_spin.setRange(1, 100000)
        self.num_designs_spin.setValue(5)
        self.num_designs_spin.setFixedWidth(96)
        design_row_layout.addWidget(self.num_designs_spin)

        design_row_layout.addWidget(QLabel("Top N designs to select"))
        self.top_n_spin = QSpinBox()
        self.top_n_spin.setRange(1, 100000)
        self.top_n_spin.setValue(1)
        self.top_n_spin.setFixedWidth(86)
        design_row_layout.addWidget(self.top_n_spin)

        design_row_layout.addWidget(QLabel("GPU"))
        self.gpu_combo = QComboBox()
        self.gpu_combo.addItems(["0", "1"])
        self.gpu_combo.setFixedWidth(70)
        design_row_layout.addWidget(self.gpu_combo)
        design_row_layout.addStretch(1)
        layout.addWidget(design_row, 4, 1, 1, 2)

        return group

    def _build_run_group(self) -> QGroupBox:
        group = QGroupBox("Parameters")
        layout = QVBoxLayout(group)
        layout.setContentsMargins(10, 12, 10, 10)
        layout.setSpacing(8)

        self.hcs_dir_label = QLabel()
        self.hcs_extracted_dir_label = QLabel()
        self._update_settings_labels()
        layout.addWidget(self.hcs_dir_label)
        layout.addWidget(self.hcs_extracted_dir_label)
        layout.addSpacing(34)

        self.generate_btn = QPushButton("Run")
        self.generate_btn.setMinimumWidth(180)
        self.generate_btn.setMinimumHeight(44)
        self.generate_btn.setCursor(Qt.PointingHandCursor)
        self.generate_btn.setStyleSheet(
            """
            QPushButton {
                background-color: #d92d20;
                color: white;
                font-size: 20px;
                font-weight: 700;
                border: 1px solid #b42318;
                border-radius: 8px;
                padding: 8px 18px;
            }
            QPushButton:hover {
                background-color: #b42318;
            }
            QPushButton:pressed {
                background-color: #912018;
            }
            QPushButton:disabled {
                background-color: #f1a7a1;
                color: #fff5f4;
                border-color: #f1a7a1;
            }
            """
        )
        self.generate_btn.clicked.connect(self._on_generate_clicked)
        layout.addWidget(self.generate_btn, 0, Qt.AlignCenter)
        layout.addStretch(1)

        return group

    def _build_input_result_group(self) -> QGroupBox:
        group = QGroupBox("Input Result")
        layout = QGridLayout(group)
        layout.setContentsMargins(8, 4, 8, 4)
        layout.setHorizontalSpacing(10)
        layout.setVerticalSpacing(3)

        self.input_result_values = {}
        columns_per_row = 3
        for index, (key, label_text) in enumerate(SCRIPT_FIELD_LABELS):
            row = index // columns_per_row
            col = index % columns_per_row
            field = QWidget()
            field_layout = QHBoxLayout(field)
            field_layout.setContentsMargins(0, 0, 0, 0)
            field_layout.setSpacing(4)

            label = QLabel(label_text)
            label.setMinimumWidth(76)
            label.setStyleSheet("color: #475467; font-weight: 700;")

            value_label = QLabel("")
            value_label.setTextInteractionFlags(Qt.TextSelectableByMouse)
            value_label.setStyleSheet("color: #101828;")
            value_label.setMinimumWidth(80)
            value_label.setWordWrap(False)

            field_layout.addWidget(label)
            field_layout.addWidget(value_label, 1)
            layout.addWidget(field, row, col)
            layout.setColumnStretch(col, 1)
            self.input_result_values[key] = value_label

        return group

    def _build_output_result_group(self) -> QGroupBox:
        group = QGroupBox("Output Result")
        layout = QVBoxLayout(group)
        layout.setContentsMargins(10, 10, 10, 10)
        layout.setSpacing(6)

        self.output_table = QTableWidget(0, len(OUTPUT_COLUMNS))
        self.output_table.setHorizontalHeaderLabels([OUTPUT_HEADER_LABELS[column] for column in OUTPUT_COLUMNS])
        self.output_table.horizontalHeader().setSectionResizeMode(QHeaderView.Stretch)
        self.output_table.verticalHeader().setVisible(False)
        self.output_table.setAlternatingRowColors(True)
        self.output_table.setSelectionBehavior(QTableWidget.SelectRows)
        self.output_table.setMinimumHeight(180)
        self.output_table.setMaximumHeight(280)
        layout.addWidget(self.output_table)

        return group

    def resizeEvent(self, event) -> None:
        super().resizeEvent(event)
        self._fit_result_tables()

    def _fit_result_tables(self) -> None:
        if hasattr(self, "output_table"):
            self.output_table.horizontalHeader().setSectionResizeMode(QHeaderView.Stretch)

    def _on_target_chain_changed(self, value: str) -> None:
        self.target_chain_other_edit.setVisible(value == "Other")

    def _get_target_chain(self) -> str:
        value = self.target_chain_combo.currentText().strip()
        if value == "Other":
            return self.target_chain_other_edit.text().strip()
        return value

    def _browse_input_pdb(self) -> None:
        file_path, _ = QFileDialog.getOpenFileName(
            self,
            "Select Input PDB",
            "",
            "Structure files (*.pdb *.cif *.mmcif);;All files (*)",
        )
        if file_path:
            self.input_pdb_edit.setText(file_path)

    def _open_settings_dialog(self) -> None:
        dialog = SettingsDialog(self.settings, self)
        if dialog.exec() == QDialog.Accepted:
            self.settings = dialog.get_settings()
            save_settings(self.settings)
            self._update_settings_labels()
            self._append_log(f"[INFO] Saved settings to {SETTINGS_PATH}")

    def _update_settings_labels(self) -> None:
        if not hasattr(self, "hcs_dir_label"):
            return
        self.hcs_dir_label.setText(f"HCS_DIR: {self.settings.get('hcs_dir', '')}")
        self.hcs_extracted_dir_label.setText(
            f"HCS_EXTRACTED_DIR: {self.settings.get('hcs_extracted_dir', '')}"
        )

    def _validate_inputs(self) -> Optional[str]:
        input_pdb = self.input_pdb_edit.text().strip()
        if not input_pdb:
            return "Please select an input PDB file."
        if not os.path.isfile(input_pdb):
            return "The selected input PDB file does not exist."
        if Path(input_pdb).suffix.lower() not in {".pdb", ".cif", ".mmcif"}:
            return "Input structure must be .pdb, .cif, or .mmcif."
        if not self._get_target_chain():
            return "Please select a target chain."
        try:
            self._parse_library_length(self.library_length_edit.text().strip())
        except ValueError as exc:
            return str(exc)
        if not self._normalize_hotspots(self.hotspot_residue_edit.text().strip(), self._get_target_chain()):
            return "Please enter at least one hotspot residue."
        if not self.settings.get("hcs_dir", "").strip():
            return "Please configure HCS_DIR in Settings."
        if not self.settings.get("hcs_extracted_dir", "").strip():
            return "Please configure HCS_EXTRACTED_DIR in Settings."
        return None

    def _parse_library_length(self, value: str) -> Tuple[int, int]:
        match = re.fullmatch(r"\s*(\d+)\s*-\s*(\d+)\s*", value)
        if not match:
            raise ValueError("Library search length must use the form 70-140.")
        minimum = int(match.group(1))
        maximum = int(match.group(2))
        if minimum <= 0 or maximum <= 0:
            raise ValueError("Library search length values must be positive.")
        if maximum < minimum:
            raise ValueError("Library search length maximum must be greater than or equal to minimum.")
        return minimum, maximum

    def _normalize_hotspots(self, raw_text: str, target_chain: str) -> List[str]:
        hotspots = []
        for token in re.split(r"[,;\s]+", raw_text.strip()):
            if not token:
                continue
            match = re.fullmatch(r"([A-Za-z]?)(\d+)", token)
            if not match:
                continue
            chain = match.group(1).upper() or target_chain
            residue_number = match.group(2)
            hotspots.append(f"{chain}{residue_number}")
        return hotspots

    def _on_generate_clicked(self) -> None:
        validation_error = self._validate_inputs()
        if validation_error:
            QMessageBox.warning(self, "Invalid Input", validation_error)
            return

        input_pdb = Path(self.input_pdb_edit.text().strip())
        target_chain = self._get_target_chain()
        min_length, max_length = self._parse_library_length(self.library_length_edit.text().strip())
        hotspots = self._normalize_hotspots(self.hotspot_residue_edit.text().strip(), target_chain)
        num_designs = int(self.num_designs_spin.value())
        top_n = int(self.top_n_spin.value())
        gpu_id = self.gpu_combo.currentText().strip()

        script_dir = input_pdb.parent / "_convex_design_gui"
        script_dir.mkdir(parents=True, exist_ok=True)
        safe_chain = re.sub(r"[^A-Za-z0-9_.-]+", "_", target_chain) or "chain"
        target_only_pdb = script_dir / f"{input_pdb.stem}_target_chain_{safe_chain}.pdb"

        try:
            convexity_metrics = calculate_paper_style_convexity(
                input_pdb=input_pdb,
                chain_id=target_chain,
                hotspot_residues=hotspots,
            )
            _write_target_chain_pdb(input_pdb, target_chain, target_only_pdb)
            target_length, target_sequence = _chain_signature(input_pdb, target_chain)
        except Exception as exc:
            QMessageBox.critical(self, "Input Processing Error", str(exc))
            return

        try:
            self.current_input_convexity = float(convexity_metrics.get("sphere convexity", 0.0) or 0.0)
        except (TypeError, ValueError):
            self.current_input_convexity = 0.0

        self.current_top_n = top_n
        self.current_selection_dir = script_dir / f"{input_pdb.stem}_5hcs_rfd" / "selected_top_designs"
        self.output_rows = []

        script_path = script_dir / f"{input_pdb.stem}_rfdiffusion_5hcs.sh"
        script_text = self._build_rfdiffusion_script(
            input_pdb=input_pdb,
            target_pdb=target_only_pdb,
            target_chain=target_chain,
            min_length=min_length,
            max_length=max_length,
            hotspots=hotspots,
            num_designs=num_designs,
            gpu_id=gpu_id,
            script_dir=script_dir,
        )
        with script_path.open("w", encoding="utf-8", newline="\n") as handle:
            handle.write(script_text)
        try:
            current_mode = script_path.stat().st_mode
            script_path.chmod(current_mode | 0o111)
        except Exception:
            pass

        self._replace_script_row(
            {
                "input file": input_pdb.name,
                "target chain": target_chain,
                "sphere convexity": self._format_metric(convexity_metrics.get("sphere convexity", "")),
                "sphere radius": self._format_metric(convexity_metrics.get("sphere radius", "")),
                "surface shape": str(convexity_metrics.get("surface shape", "")),
                "library length": f"{min_length}-{max_length}",
                "hotspot residues": ",".join(hotspots),
                "number of designs": str(num_designs),
                "gpu": gpu_id,
            }
        )
        self._replace_output_row(
            {
                "pdb file": "(pending)",
                "sphere convexity": "",
                "length": "",
                "status": "Completed",
            }
        )
        self.log_edit.setPlainText(script_text)
        self.progress_bar.setValue(0)
        self._append_log(f"\n[INFO] Wrote RFdiffusion script to {script_path}")
        self._append_log(f"[INFO] Wrote target-only PDB to {target_only_pdb}")
        self._start_script_run(
            script_path,
            script_dir / f"{input_pdb.stem}_5hcs_rfd" / "rfd_result",
            num_designs,
            target_chain,
            target_sequence,
            target_length,
        )

    def _format_metric(self, value: object) -> str:
        if value in ("", None):
            return ""
        try:
            return f"{float(value):.4f}"
        except (TypeError, ValueError):
            return str(value)

    def _build_rfdiffusion_script(
        self,
        input_pdb: Path,
        target_pdb: Path,
        target_chain: str,
        min_length: int,
        max_length: int,
        hotspots: List[str],
        num_designs: int,
        gpu_id: str,
        script_dir: Path,
    ) -> str:
        hcs_dir = self.settings["hcs_dir"]
        hcs_extracted_dir = self.settings["hcs_extracted_dir"]
        run_root = script_dir / f"{input_pdb.stem}_5hcs_rfd"
        hotspot_text = "[" + ",".join(hotspots) + "]"
        return f"""#!/bin/bash
set -euo pipefail

module load miniforge3/
source "$(conda info --base)/etc/profile.d/conda.sh"
set +u
conda deactivate
conda activate SE3nv
set -u

export CUDA_VISIBLE_DEVICES={gpu_id}

RFDIFFUSION_DIR="${{RFDIFFUSION_DIR:-/home/cs/RFdiffusion}}"
TARGET_PDB="{target_pdb.as_posix()}"
TARGET_CHAIN="{target_chain}"
FIVE_HCS_DIR="{Path(hcs_dir).as_posix()}"
FIVE_HCS_EXTRACTED_DIR="{Path(hcs_extracted_dir).as_posix()}"
RUN_ROOT="{run_root.as_posix()}"

TARGET_NAME="{input_pdb.stem}"
TARGET_HOTSPOTS="{hotspot_text}"
NUM_DESIGNS={num_designs}
MIN_SCAFFOLD_RESIDUES={min_length}
MAX_SCAFFOLD_RESIDUES={max_length}

SCAFFOLD_PDB_DIR="${{RUN_ROOT}}/5hcs_pdbs"
SCAFFOLD_GUIDE_DIR="${{RUN_ROOT}}/5hcs_scaffoldguided"
TARGET_GUIDE_DIR="${{RUN_ROOT}}/target_scaffoldguided"
OUTPUT_PREFIX="${{RUN_ROOT}}/rfd_result/${{TARGET_NAME}}_5hcs"

rm -rf "${{SCAFFOLD_PDB_DIR}}" "${{SCAFFOLD_GUIDE_DIR}}"
mkdir -p "${{SCAFFOLD_PDB_DIR}}" "${{SCAFFOLD_GUIDE_DIR}}" "${{TARGET_GUIDE_DIR}}" "$(dirname "${{OUTPUT_PREFIX}}")"

if [[ ! -d "${{FIVE_HCS_DIR}}" ]]; then
  echo "5HCS directory does not exist: ${{FIVE_HCS_DIR}}" >&2
  exit 1
fi

if [[ ! -d "${{FIVE_HCS_EXTRACTED_DIR}}" ]]; then
  echo "Extracted 5HCS directory does not exist: ${{FIVE_HCS_EXTRACTED_DIR}}" >&2
  exit 1
fi

set +o pipefail
find "${{FIVE_HCS_EXTRACTED_DIR}}" -maxdepth 4 -type f | head -30
set -o pipefail

CACHE_PARENT="${{FIVE_HCS_EXTRACTED_DIR}}/.convex_design_gui_cache"
if mkdir -p "${{CACHE_PARENT}}" 2>/dev/null; then
  CACHE_DIR="${{CACHE_PARENT}}"
else
  CACHE_KEY="$(printf '%s' "${{FIVE_HCS_EXTRACTED_DIR}}" | cksum | awk '{{print $1}}')"
  CACHE_DIR="${{HOME}}/.cache/convex_design_gui/5hcs_${{CACHE_KEY}}"
  mkdir -p "${{CACHE_DIR}}"
fi
SCAFFOLD_INDEX="${{CACHE_DIR}}/5hcs_scaffold_lengths.tsv"

if [[ "${{REBUILD_HCS_INDEX:-0}}" == "1" ]]; then
  rm -f "${{SCAFFOLD_INDEX}}"
fi

if [[ ! -s "${{SCAFFOLD_INDEX}}" ]]; then
  echo "Building 5HCS scaffold length index: ${{SCAFFOLD_INDEX}}"
  index_tmp="${{SCAFFOLD_INDEX}}.tmp.$$"
  indexed_count=0
  while IFS= read -r -d '' scaffold_path; do
    lower_path="$(echo "${{scaffold_path}}" | tr '[:upper:]' '[:lower:]')"
    if [[ "${{lower_path}}" == *.gz ]]; then
      residue_count="$(
        gzip -cd "${{scaffold_path}}" | awk '
          $1 == "ATOM" && $3 == "CA" {{
            key = $5 ":" $6
            seen[key] = 1
          }}
          END {{
            count = 0
            for (key in seen) {{
              count++
            }}
            print count
          }}
        '
      )"
    else
      residue_count="$(
        awk '
          $1 == "ATOM" && $3 == "CA" {{
            key = $5 ":" $6
            seen[key] = 1
          }}
          END {{
            count = 0
            for (key in seen) {{
              count++
            }}
            print count
          }}
        ' "${{scaffold_path}}"
      )"
    fi
    scaffold_name="${{scaffold_path#${{FIVE_HCS_EXTRACTED_DIR}}/}}"
    scaffold_name="$(echo "${{scaffold_name}}" | tr '/ ' '__')"
    scaffold_name="${{scaffold_name%.gz}}"
    scaffold_name="${{scaffold_name%.*}}.pdb"
    printf '%s\t%s\t%s\n' "${{residue_count}}" "${{scaffold_name}}" "${{scaffold_path}}" >> "${{index_tmp}}"
    indexed_count=$((indexed_count + 1))
    if (( indexed_count % 1000 == 0 )); then
      echo "Indexed ${{indexed_count}} scaffolds."
    fi
  done < <(
    find "${{FIVE_HCS_EXTRACTED_DIR}}" \\
      -type f \\
      \\( -iname "*.pdb" -o -iname "*.pdb.gz" \\) \\
      -print0
  )
  mv "${{index_tmp}}" "${{SCAFFOLD_INDEX}}"
  echo "Indexed ${{indexed_count}} scaffold structure files."
else
  echo "Using cached 5HCS scaffold length index: ${{SCAFFOLD_INDEX}}"
fi

selected_count=0
rejected_count=0
processed_count=0

while IFS=$'\t' read -r residue_count scaffold_name scaffold_path; do
  if [[ -z "${{residue_count}}" || -z "${{scaffold_name}}" || -z "${{scaffold_path}}" ]]; then
    continue
  fi
  processed_count=$((processed_count + 1))
  if [[ "${{residue_count}}" -lt "${{MIN_SCAFFOLD_RESIDUES}}" || "${{residue_count}}" -gt "${{MAX_SCAFFOLD_RESIDUES}}" ]]; then
    rejected_count=$((rejected_count + 1))
    continue
  fi

  scaffold_out="${{SCAFFOLD_PDB_DIR}}/${{scaffold_name}}"
  if [[ -e "${{scaffold_out}}" || -L "${{scaffold_out}}" ]]; then
    continue
  fi
  lower_path="$(echo "${{scaffold_path}}" | tr '[:upper:]' '[:lower:]')"
  if [[ "${{lower_path}}" == *.gz ]]; then
    gzip -cd "${{scaffold_path}}" > "${{scaffold_out}}"
  else
    ln -s "${{scaffold_path}}" "${{scaffold_out}}"
  fi
  selected_count=$((selected_count + 1))
  if (( processed_count % 1000 == 0 )); then
    echo "Processed ${{processed_count}} scaffolds. Selected ${{selected_count}}, rejected ${{rejected_count}}."
  fi
done < "${{SCAFFOLD_INDEX}}"

echo "Selected ${{selected_count}} scaffold structure files."
if ! find -L "${{SCAFFOLD_PDB_DIR}}" -type f -name "*.pdb" -print -quit | grep -q .; then
  echo "No 5HCS scaffold structure files passed the current scan and length filter." >&2
  exit 1
fi

python "${{RFDIFFUSION_DIR}}/helper_scripts/make_secstruc_adj.py" \\
  --pdb_dir "${{SCAFFOLD_PDB_DIR}}" \\
  --out_dir "${{SCAFFOLD_GUIDE_DIR}}"

python "${{RFDIFFUSION_DIR}}/helper_scripts/make_secstruc_adj.py" \\
  --input_pdb "${{TARGET_PDB}}" \\
  --out_dir "${{TARGET_GUIDE_DIR}}"

TARGET_STEM="$(basename "${{TARGET_PDB}}")"
TARGET_STEM="${{TARGET_STEM%.*}}"
TARGET_SS="${{TARGET_GUIDE_DIR}}/${{TARGET_STEM}}_ss.pt"
TARGET_ADJ="${{TARGET_GUIDE_DIR}}/${{TARGET_STEM}}_adj.pt"

set +e
python "${{RFDIFFUSION_DIR}}/scripts/run_inference.py" \\
  scaffoldguided.target_path="${{TARGET_PDB}}" \\
  inference.output_prefix="${{OUTPUT_PREFIX}}" \\
  scaffoldguided.scaffoldguided=True \\
  "ppi.hotspot_res=${{TARGET_HOTSPOTS}}" \\
  scaffoldguided.target_pdb=True \\
  scaffoldguided.target_ss="${{TARGET_SS}}" \\
  scaffoldguided.target_adj="${{TARGET_ADJ}}" \\
  scaffoldguided.scaffold_dir="${{SCAFFOLD_GUIDE_DIR}}" \\
  scaffoldguided.mask_loops=False \\
  inference.num_designs="${{NUM_DESIGNS}}" \\
  denoiser.noise_scale_ca=0.5 \\
  denoiser.noise_scale_frame=0.5
rf_exit_code=$?
set -e

output_pdb_count="$(find "$(dirname "${{OUTPUT_PREFIX}}")" -type f -name "*.pdb" -not -path "*/traj/*" | wc -l)"
echo "Output PDB count: ${{output_pdb_count}}"

if [[ "${{rf_exit_code}}" -ne 0 ]]; then
  if [[ "${{rf_exit_code}}" -eq 139 && "${{output_pdb_count}}" -gt 0 ]]; then
    echo "RFdiffusion exited with segmentation fault after writing output PDB files."
    echo "Treating this run as completed with warning."
    exit 0
  fi
  echo "RFdiffusion failed with exit code ${{rf_exit_code}}." >&2
  exit "${{rf_exit_code}}"
fi
"""

    def _replace_script_row(self, row: Dict[str, str]) -> None:
        for key, value_label in self.input_result_values.items():
            value = row.get(key, "")
            value_label.setText(value)
            value_label.setToolTip(value)

    def _replace_output_row(self, row: Dict[str, str]) -> None:
        self.output_table.setRowCount(0)
        row_index = self.output_table.rowCount()
        self.output_table.insertRow(row_index)
        for column_index, column_name in enumerate(OUTPUT_COLUMNS):
            value = row.get(column_name, "")
            item = QTableWidgetItem(value)
            item.setToolTip(value)
            if column_name in {"sphere convexity", "convexity difference", "length"}:
                item.setTextAlignment(Qt.AlignRight | Qt.AlignVCenter)
            self.output_table.setItem(row_index, column_index, item)

    def _set_output_status(self, status: str) -> None:
        if self.output_rows:
            for row in self.output_rows:
                row["status"] = status
            self._sort_output_table_by_difference()
            return

        if self.output_table.rowCount() == 0:
            return

        status_index = OUTPUT_COLUMNS.index("status")
        for row_index in range(self.output_table.rowCount()):
            item = self.output_table.item(row_index, status_index)
            if item is None:
                item = QTableWidgetItem()
                self.output_table.setItem(row_index, status_index, item)
            item.setText(status)
            item.setToolTip(status)

    def _start_script_run(
        self,
        script_path: Path,
        output_dir: Path,
        expected_designs: int,
        target_chain: str,
        target_sequence: str,
        target_length: int,
    ) -> None:
        self.generate_btn.setEnabled(False)
        self._set_output_status("running")
        self.run_thread = QThread(self)
        self.run_worker = ScriptRunWorker(
            script_path,
            output_dir,
            expected_designs,
            target_chain,
            target_sequence,
            target_length,
        )
        self.run_worker.moveToThread(self.run_thread)

        self.run_thread.started.connect(self.run_worker.run)
        self.run_worker.signals.log.connect(self._append_log)
        self.run_worker.signals.error.connect(self._on_run_error)
        self.run_worker.signals.status.connect(self._set_output_status)
        self.run_worker.signals.output_count.connect(self._set_output_pdb_count)
        self.run_worker.signals.output_pdb.connect(self._add_output_pdb_row)
        self.run_worker.signals.progress.connect(self.progress_bar.setValue)
        self.run_worker.signals.finished.connect(self._on_run_finished)
        self.run_worker.signals.finished.connect(self.run_thread.quit)
        self.run_worker.signals.finished.connect(self.run_worker.deleteLater)
        self.run_thread.finished.connect(self.run_thread.deleteLater)
        self.run_thread.start()

    def _on_run_error(self, error_text: str) -> None:
        self._append_log("[ERROR] " + error_text)
        QMessageBox.critical(self, "Error", error_text)

    def _on_run_finished(self, return_code: int) -> None:
        self.generate_btn.setEnabled(True)
        if return_code == 0:
            self._copy_top_designs()
        if return_code == 0:
            self._append_log("[INFO] RFdiffusion job finished.")
        elif return_code == 139:
            self._append_log("[ERROR] RFdiffusion exited with segmentation fault after process shutdown.")
        else:
            self._append_log(f"[ERROR] RFdiffusion job exited with code {return_code}.")
        self.run_worker = None
        self.run_thread = None

    def _set_output_pdb_count(self, count: int) -> None:
        self.progress_bar.setValue(max(self.progress_bar.value(), min(100, int(count))))

    def _add_output_pdb_row(self, row: Dict[str, object]) -> None:
        output_convexity = row.get("sphere convexity", "")
        try:
            output_convexity_value = float(output_convexity)
            difference = abs(abs(self.current_input_convexity) - abs(output_convexity_value))
            shape_match = "Yes" if self.current_input_convexity * output_convexity_value < 0.0 else "No"
        except (TypeError, ValueError):
            difference = ""
            shape_match = ""
        normalized_row = dict(row)
        normalized_row["convexity difference"] = difference
        normalized_row["shape match"] = shape_match
        row_path = str(normalized_row.get("path", ""))
        existing_row = None
        if row_path:
            for current_row in self.output_rows:
                if str(current_row.get("path", "")) == row_path:
                    existing_row = current_row
                    break

        if existing_row is None:
            self.output_rows.append(normalized_row)
        else:
            existing_row.update(normalized_row)
        self._sort_output_table_by_difference()

    def _sort_output_table_by_difference(self) -> None:
        sorted_rows = sorted(
            self.output_rows,
            key=self._output_rank_key,
        )
        self._render_output_rows(sorted_rows)

    def _render_output_rows(self, rows: List[Dict[str, object]]) -> None:
        self.output_table.setRowCount(0)
        for row in rows:
            row_index = self.output_table.rowCount()
            self.output_table.insertRow(row_index)
            for column_index, column_name in enumerate(OUTPUT_COLUMNS):
                value = row.get(column_name, "")
                if column_name in {"sphere convexity", "convexity difference"}:
                    value = self._format_metric(value)
                item = QTableWidgetItem(str(value))
                item.setToolTip(str(value))
                if column_name in {"sphere convexity", "convexity difference", "length"}:
                    item.setTextAlignment(Qt.AlignRight | Qt.AlignVCenter)
                self.output_table.setItem(row_index, column_index, item)

    def _copy_top_designs(self) -> None:
        if self.current_selection_dir is None or not self.output_rows:
            return
        sorted_rows = sorted(
            self.output_rows,
            key=self._output_rank_key,
        )
        selected_rows = sorted_rows[: self.current_top_n]
        selected_paths = {str(row.get("path", "")) for row in selected_rows}
        for row in self.output_rows:
            if str(row.get("path", "")) in selected_paths:
                row["status"] = "Selected"
            else:
                row["status"] = row.get("status", "") or "Completed"
        self.current_selection_dir.mkdir(parents=True, exist_ok=True)
        for row in selected_rows:
            source_path = Path(str(row.get("path", "")))
            if not source_path.is_file():
                continue
            target_path = self.current_selection_dir / source_path.name
            shutil.copy2(source_path, target_path)
        summary_path = self.current_selection_dir / "selected_top_designs.csv"
        with summary_path.open("w", encoding="utf-8", newline="") as handle:
            writer = csv.DictWriter(
                handle,
                fieldnames=[
                    "rank",
                    "pdb file",
                    "source path",
                    "sphere convexity",
                    "convexity difference",
                    "shape match",
                    "length",
                    "status",
                ],
            )
            writer.writeheader()
            for rank, row in enumerate(selected_rows, start=1):
                writer.writerow(
                    {
                        "rank": rank,
                        "pdb file": row.get("pdb file", ""),
                        "source path": row.get("path", ""),
                        "sphere convexity": row.get("sphere convexity", ""),
                        "convexity difference": row.get("convexity difference", ""),
                        "shape match": row.get("shape match", ""),
                        "length": row.get("length", ""),
                        "status": "Selected",
                    }
                )
        self._sort_output_table_by_difference()
        self._append_log(
            f"[INFO] Copied {len(selected_rows)} top designs to {self.current_selection_dir}"
        )
        self._append_log(f"[INFO] Wrote top design summary to {summary_path}")

    def _output_rank_key(self, row: Dict[str, object]) -> Tuple[int, float]:
        shape_match = 0 if str(row.get("shape match", "")) == "Yes" else 1
        try:
            difference = float(row.get("convexity difference", 999999.0) or 999999.0)
        except (TypeError, ValueError):
            difference = 999999.0
        return shape_match, difference

    def _append_log(self, message: str) -> None:
        self.log_edit.appendPlainText(str(message))
        output_match = re.search(r"Output PDB count:\s*(\d+)", str(message))
        if output_match and self.output_table.rowCount() > 0:
            self._set_output_pdb_count(int(output_match.group(1)))
        cursor = self.log_edit.textCursor()
        cursor.movePosition(cursor.MoveOperation.End)
        self.log_edit.setTextCursor(cursor)

    def closeEvent(self, event) -> None:
        if self.run_worker is not None and self.generate_btn is not None and not self.generate_btn.isEnabled():
            response = QMessageBox.question(
                self,
                "Exit",
                "A job is still running. Close the GUI and terminate active jobs?",
                QMessageBox.Yes | QMessageBox.No,
                QMessageBox.No,
            )
            if response != QMessageBox.Yes:
                event.ignore()
                return
            self.run_worker.request_cancel()
        event.accept()


def _ensure_writable_runtime_dir() -> None:
    if os.name != "posix":
        return
    current_runtime_dir = os.environ.get("XDG_RUNTIME_DIR", "").strip()
    if current_runtime_dir:
        current_path = Path(current_runtime_dir)
        if current_path.is_dir() and os.access(str(current_path), os.W_OK):
            return

    user_name = os.environ.get("USER", "user")
    fallback_path = Path("/tmp") / f"runtime-{user_name}-convex-design-gui"
    try:
        fallback_path.mkdir(mode=0o700, parents=True, exist_ok=True)
        fallback_path.chmod(0o700)
        os.environ["XDG_RUNTIME_DIR"] = str(fallback_path)
    except Exception:
        pass


def run() -> None:
    _ensure_writable_runtime_dir()
    app = QApplication(sys.argv)
    window = MainWindow()
    window.show()
    sys.exit(app.exec())
