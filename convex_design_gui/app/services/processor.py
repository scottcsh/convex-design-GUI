import csv
import math
import random
from dataclasses import dataclass
from pathlib import Path
from typing import Callable, Dict, Iterable, List, Optional, Sequence, Set, Tuple

from Bio.PDB import MMCIFParser, PDBIO, PDBParser, Select, is_aa


SUPPORTED_EXTENSIONS = {".pdb", ".cif", ".mmcif"}


@dataclass(frozen=True)
class ResiduePoint:
    chain_id: str
    residue_number: int
    residue_name: str
    coord: Tuple[float, float, float]


class SelectedResidues(Select):
    def __init__(self, selected_keys: Set[Tuple[str, int]]) -> None:
        super().__init__()
        self.selected_keys = selected_keys

    def accept_residue(self, residue) -> bool:
        parent = residue.get_parent()
        key = (parent.id, int(residue.id[1]))
        return key in self.selected_keys


def process_directory(
    input_dir: str,
    output_csv: str,
    chain_id: str,
    residue_ranges: str = "",
    exposure_percentile: float = 70.0,
    convexity_cutoff: float = 0.0,
    neighbor_count: int = 8,
    write_partial_pdb: bool = True,
    log_callback: Optional[Callable[[str], None]] = None,
    result_row_callback: Optional[Callable[[Dict[str, object]], None]] = None,
    total_progress_callback: Optional[Callable[[int], None]] = None,
    current_progress_callback: Optional[Callable[[int], None]] = None,
    status_callback: Optional[Callable[[str], None]] = None,
    should_cancel_callback: Optional[Callable[[], bool]] = None,
) -> List[Dict[str, object]]:
    input_path = Path(input_dir)
    output_path = Path(output_csv)
    structure_files = _find_structure_files(input_path)

    if not structure_files:
        raise ValueError("No .pdb, .cif, or .mmcif files were found in the selected input directory.")
    if not chain_id.strip():
        raise ValueError("A chain must be selected.")

    output_path.parent.mkdir(parents=True, exist_ok=True)
    partial_dir = output_path.parent / "_convex_design_gui" / "partials"
    if write_partial_pdb:
        partial_dir.mkdir(parents=True, exist_ok=True)

    selected_ranges = parse_residue_ranges(residue_ranges, default_chain=chain_id.strip())
    rows: List[Dict[str, object]] = []

    _emit(total_progress_callback, 0)
    for file_index, structure_path in enumerate(structure_files, start=1):
        _raise_if_cancelled(should_cancel_callback)
        _emit_status(status_callback, f"Processing {structure_path.name}")
        _emit(current_progress_callback, 0)
        _log(log_callback, f"[INFO] Processing {structure_path.name}")

        row = process_structure(
            structure_path=structure_path,
            chain_id=chain_id.strip(),
            selected_ranges=selected_ranges,
            exposure_percentile=exposure_percentile,
            convexity_cutoff=convexity_cutoff,
            neighbor_count=neighbor_count,
            partial_dir=partial_dir if write_partial_pdb else None,
        )
        rows.append(row)
        _emit_result(result_row_callback, row)

        progress = int(file_index / max(len(structure_files), 1) * 100)
        _emit(current_progress_callback, 100)
        _emit(total_progress_callback, progress)

    write_results_csv(output_path, rows)
    _log(log_callback, f"[INFO] Wrote {len(rows)} rows to {output_path}")
    _emit_status(status_callback, "Finished")
    _emit(total_progress_callback, 100)
    _emit(current_progress_callback, 100)
    return rows


def process_structure(
    structure_path: Path,
    chain_id: str,
    selected_ranges: Set[Tuple[str, int]],
    exposure_percentile: float,
    convexity_cutoff: float,
    neighbor_count: int,
    partial_dir: Optional[Path],
) -> Dict[str, object]:
    structure = _load_structure(structure_path)
    model = structure[0]
    if chain_id not in model:
        available = ", ".join(chain.id for chain in model)
        raise ValueError(f"Chain '{chain_id}' was not found in {structure_path.name}. Available chains: {available or '(none)'}")

    residue_points = _collect_residue_points(model[chain_id])
    if not residue_points:
        raise ValueError(f"Chain '{chain_id}' in {structure_path.name} has no protein residues with CA atoms.")

    ranged_points = [
        point for point in residue_points if not selected_ranges or (point.chain_id, point.residue_number) in selected_ranges
    ]
    if not ranged_points:
        raise ValueError(f"No residues matched the selected range for {structure_path.name}.")

    centroid = _centroid([point.coord for point in residue_points])
    radial_distances = {point.residue_number: _distance(point.coord, centroid) for point in residue_points}
    exposure_threshold = _percentile(list(radial_distances.values()), exposure_percentile)
    convexity_scores = _convexity_scores(residue_points, radial_distances, neighbor_count)

    selected_points = []
    for point in ranged_points:
        radial_distance = radial_distances[point.residue_number]
        convexity = convexity_scores[point.residue_number]
        if radial_distance >= exposure_threshold and convexity >= convexity_cutoff:
            selected_points.append(point)

    selected_keys = {(point.chain_id, point.residue_number) for point in selected_points}
    partial_path = ""
    if partial_dir is not None and selected_keys:
        partial_path = str(partial_dir / f"{structure_path.stem}_convex_partial.pdb")
        io = PDBIO()
        io.set_structure(structure)
        io.save(partial_path, SelectedResidues(selected_keys))

    selected_residue_text = ",".join(f"{point.chain_id}{point.residue_number}" for point in selected_points)
    average_convexity = _average([convexity_scores[point.residue_number] for point in selected_points])
    average_exposure = _average([radial_distances[point.residue_number] for point in selected_points])

    return {
        "input file name": structure_path.name,
        "chain": chain_id,
        "candidate residues": len(ranged_points),
        "convex residues": len(selected_points),
        "mean convexity": average_convexity,
        "mean radial exposure": average_exposure,
        **_calculate_sphere_surface_metrics(
            model=model,
            chain_id=chain_id,
            selected_ranges=selected_ranges,
            fallback_points=ranged_points,
        ),
        "selected residues": selected_residue_text,
        "partial pdb": partial_path,
    }


def parse_residue_ranges(raw_text: str, default_chain: str) -> Set[Tuple[str, int]]:
    selected: Set[Tuple[str, int]] = set()
    text = raw_text.strip()
    if not text:
        return selected

    for raw_token in text.replace(";", ",").split(","):
        token = raw_token.strip()
        if not token:
            continue
        chain = default_chain
        residue_part = token
        if token[0].isalpha():
            chain = token[0]
            residue_part = token[1:]
        if "-" in residue_part:
            start_text, end_text = residue_part.split("-", 1)
            start = int(start_text.strip())
            end = int(end_text.strip())
            if end < start:
                start, end = end, start
            for residue_number in range(start, end + 1):
                selected.add((chain, residue_number))
        else:
            selected.add((chain, int(residue_part)))
    return selected


def write_results_csv(output_path: Path, rows: List[Dict[str, object]]) -> None:
    columns = [
        "input file name",
        "chain",
        "candidate residues",
        "convex residues",
        "mean convexity",
        "mean radial exposure",
        "sphere radius",
        "sphere convexity",
        "surface shape",
        "sphere inliers",
        "selected residues",
        "partial pdb",
    ]
    with output_path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=columns)
        writer.writeheader()
        for row in rows:
            writer.writerow({column: _format_value(row.get(column, "")) for column in columns})


def _find_structure_files(input_path: Path) -> List[Path]:
    return sorted(path for path in input_path.iterdir() if path.is_file() and path.suffix.lower() in SUPPORTED_EXTENSIONS)


def _load_structure(source_path: Path):
    suffix = source_path.suffix.lower()
    if suffix == ".pdb":
        parser = PDBParser(QUIET=True)
    elif suffix in {".cif", ".mmcif"}:
        parser = MMCIFParser(QUIET=True)
    else:
        raise ValueError(f"Unsupported structure format: {source_path.suffix}")
    return parser.get_structure(source_path.stem, str(source_path))


def _collect_residue_points(chain) -> List[ResiduePoint]:
    points: List[ResiduePoint] = []
    for residue in chain:
        if residue.id[0] != " " or not is_aa(residue, standard=False) or "CA" not in residue:
            continue
        atom = residue["CA"]
        coord = atom.get_coord()
        points.append(
            ResiduePoint(
                chain_id=chain.id,
                residue_number=int(residue.id[1]),
                residue_name=residue.get_resname(),
                coord=(float(coord[0]), float(coord[1]), float(coord[2])),
            )
        )
    return points


def _calculate_sphere_surface_metrics(
    model,
    chain_id: str,
    selected_ranges: Set[Tuple[str, int]],
    fallback_points: Sequence[ResiduePoint],
) -> Dict[str, object]:
    chain = model[chain_id]
    chain_heavy_coords = _collect_heavy_atom_coords(chain, selected_ranges=set())
    patch_heavy_coords = _collect_heavy_atom_coords(chain, selected_ranges=selected_ranges)
    if not patch_heavy_coords:
        patch_heavy_coords = _collect_heavy_atom_coords(
            chain,
            selected_ranges={(point.chain_id, point.residue_number) for point in fallback_points},
        )

    if len(patch_heavy_coords) < 4 or len(chain_heavy_coords) < 4:
        return {
            "sphere radius": "",
            "sphere convexity": "",
            "surface shape": "",
            "sphere inliers": "",
        }

    sphere = _ransac_sphere_fit(
        coords=patch_heavy_coords,
        threshold=1.0,
        max_iterations=2000,
        seed=17,
    )
    if sphere is None:
        return {
            "sphere radius": "",
            "sphere convexity": "",
            "surface shape": "",
            "sphere inliers": "",
        }

    sphere_center, radius, inlier_count = sphere
    chain_center = _centroid(chain_heavy_coords)
    patch_center = _centroid(patch_heavy_coords)
    protein_vector = _vector(patch_center, chain_center)
    fitted_vector = _vector(patch_center, sphere_center)
    direction_score = _dot(protein_vector, fitted_vector)
    surface_shape = "convex" if direction_score >= 0.0 else "concave"
    magnitude = 1.0 / radius if radius > 0.0 else ""
    signed_convexity = -magnitude if direction_score >= 0.0 else magnitude
    return {
        "sphere radius": radius,
        "sphere convexity": signed_convexity,
        "surface shape": surface_shape,
        "sphere inliers": inlier_count,
    }


def _collect_heavy_atom_coords(chain, selected_ranges: Set[Tuple[str, int]]) -> List[Tuple[float, float, float]]:
    coords: List[Tuple[float, float, float]] = []
    for residue in chain:
        if residue.id[0] != " " or not is_aa(residue, standard=False):
            continue
        key = (chain.id, int(residue.id[1]))
        if selected_ranges and key not in selected_ranges:
            continue
        for atom in residue:
            element = str(getattr(atom, "element", "")).strip().upper()
            if element == "H":
                continue
            coord = atom.get_coord()
            coords.append((float(coord[0]), float(coord[1]), float(coord[2])))
    return coords


def _ransac_sphere_fit(
    coords: Sequence[Tuple[float, float, float]],
    threshold: float,
    max_iterations: int,
    seed: int,
) -> Optional[Tuple[Tuple[float, float, float], float, int]]:
    if len(coords) < 4:
        return None

    rng = random.Random(seed)
    best_center: Optional[Tuple[float, float, float]] = None
    best_radius = 0.0
    best_inliers = -1
    iterations = min(max_iterations, max(200, len(coords) * 20))

    for _ in range(iterations):
        sample = rng.sample(list(coords), 4)
        sphere = _sphere_from_four_points(sample)
        if sphere is None:
            continue
        center, radius = sphere
        if radius <= 0.0 or not math.isfinite(radius):
            continue
        inliers = sum(1 for coord in coords if abs(_distance(coord, center) - radius) <= threshold)
        if inliers > best_inliers:
            best_center = center
            best_radius = radius
            best_inliers = inliers

    if best_center is None:
        return None
    return best_center, best_radius, best_inliers


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
    radius_squared = (
        center[0] * center[0]
        + center[1] * center[1]
        + center[2] * center[2]
        - d_value
    )
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


def _convexity_scores(
    residue_points: Sequence[ResiduePoint],
    radial_distances: Dict[int, float],
    neighbor_count: int,
) -> Dict[int, float]:
    scores: Dict[int, float] = {}
    effective_neighbor_count = max(1, min(neighbor_count, max(len(residue_points) - 1, 1)))
    for point in residue_points:
        neighbors = sorted(
            (other for other in residue_points if other.residue_number != point.residue_number),
            key=lambda other: _distance(point.coord, other.coord),
        )[:effective_neighbor_count]
        neighbor_average = _average(radial_distances[other.residue_number] for other in neighbors)
        scores[point.residue_number] = radial_distances[point.residue_number] - neighbor_average
    return scores


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


def _percentile(values: Sequence[float], percentile: float) -> float:
    if not values:
        return 0.0
    ordered = sorted(values)
    clamped = max(0.0, min(100.0, float(percentile)))
    index = (len(ordered) - 1) * clamped / 100.0
    lower = int(math.floor(index))
    upper = int(math.ceil(index))
    if lower == upper:
        return ordered[lower]
    fraction = index - lower
    return ordered[lower] * (1.0 - fraction) + ordered[upper] * fraction


def _average(values: Iterable[float]) -> float:
    collected = list(values)
    if not collected:
        return 0.0
    return sum(collected) / len(collected)


def _format_value(value: object) -> object:
    if isinstance(value, float):
        return f"{value:.4f}"
    return value


def _emit(callback: Optional[Callable[[int], None]], value: int) -> None:
    if callback is not None:
        callback(max(0, min(100, int(value))))


def _emit_status(callback: Optional[Callable[[str], None]], value: str) -> None:
    if callback is not None:
        callback(str(value))


def _emit_result(callback: Optional[Callable[[Dict[str, object]], None]], row: Dict[str, object]) -> None:
    if callback is not None:
        callback(dict(row))


def _log(callback: Optional[Callable[[str], None]], message: str) -> None:
    if callback is not None:
        callback(str(message))


def _raise_if_cancelled(callback: Optional[Callable[[], bool]]) -> None:
    if callback is not None and callback():
        raise RuntimeError("Cancelled")
