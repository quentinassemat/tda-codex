import shlex
import numpy as np

class CodexSample:
    """Class representing a single CODEX sample."""
    cells_by_phenotype : dict[str,list[tuple[int, int]]] #TODO: maybe create a Cell type and store a list of these instead?
    def __init__(self) -> None:
        self.cells_by_phenotype = dict()

    def add_cell(self, cell_x: int, cell_y: int, phenotype: str):
        """Adds a cell with the given phenotype to the sample.

        Args:
            cell_x (int): X coordinate of the cell to add
            cell_y (int): Y coordinate of the cell to add
            phenotype (str): _description_
        """
        if phenotype not in self.cells_by_phenotype:
            self.cells_by_phenotype[phenotype] = []
        self.cells_by_phenotype[phenotype].append((cell_x, cell_y))

    def get_cells(self, phenotype: str):
        """Returns a list of cells of the given phenotype in this sample, or an empty list if there are none.

        Args:
            phenotype (str): The phenotype to return.

        Returns:
            tuple[list[int], list[int]]: The list of x- and y-coordinates of the cells of the given phenotype in this sample. Empty if no such cells exist.
        """
        if phenotype in self.cells_by_phenotype:
            cells = self.cells_by_phenotype[phenotype]
            xs = [pair[0] for pair in cells]
            ys = [pair[1] for pair in cells]
            return xs,ys
        else:
            return [], []

    def get_subsample(self, x: int, y: int, width: int, height: int):
        """Returns a rectangular subsample of this sample.

        Args:
            x (int): The x-coordinate of the upper-left corner of the subsample
            y (int): The y-coordinate of the upper-left corner of the subsample
            width (int): The width of the subsample
            height (int): The height of the subsample

        Returns:
            CodexSample: A sample consisting only of the cells in the region specified.
        """
        subsample = CodexSample()
        for phenotype, cell_list in self.cells_by_phenotype.items():
            for cell_x, cell_y in cell_list:
                if cell_x >= x and cell_x < x+width and cell_y >= y and cell_y < y+height:
                    subsample.add_cell(cell_x, cell_y, phenotype)
        return subsample

    def get_num_cells(self, phenotype: str):
        return len(self.get_cells(phenotype))
        
    def get_top_phenotypes(self):
        return sorted(self.cells_by_phenotype.keys(), key=self.get_num_cells, reverse=True)

class CodexData:
    """Class providing access to a set of CODEX samples."""
    known_phenotypes : list[str]
    samples : dict[str, CodexSample]
    def __init__(self) -> None:
        self.samples = dict()
        self.known_phenotypes = []

    def add_cell(self, sample, xx, yy, phenotype):
        """Adds a cell to a sample. Creates a new sample if one does not already exist with the given name.

        Args:
            sample (_type_): _description_
            xx (_type_): _description_
            yy (_type_): _description_
            phenotype (_type_): _description_
        """
        if phenotype not in self.known_phenotypes:
            self.known_phenotypes.append(phenotype)
        if sample not in self.samples.keys():
            self.samples[sample] = CodexSample()
        self.samples[sample].add_cell(xx, yy, phenotype)

    def get_num_cells(self, phenotype: str):
        total = 0
        for sample in self.samples.values():
            total += len(sample.get_cells(phenotype))
        return total

    def get_top_phenotypes(self):
        return sorted(self.known_phenotypes, key=self.get_num_cells, reverse=True)


def read_cluster_annotations(fname) -> dict[int, str]:
    cluster_annotations = dict()
    with open(fname) as f:
        for line in f:
            key, val = line.split(",")
            if key in cluster_annotations.keys():
                raise ValueError(f"Duplicate cluster annotation '{key}'")
            cluster_annotations[key] = val.strip()
    return cluster_annotations

def read_raw_spleen_csv(raw_spleen_csv_path: str, spleen_annotations_path: str) -> CodexData:
    cluster_annotations = read_cluster_annotations(spleen_annotations_path)
    x_tile_size = 1344 #TODO: maybe make this not be hardcoded?
    y_tile_size = 1008

    with open(raw_spleen_csv_path) as f:
        codex = CodexData()
        headers = f.readline().split(",")
        for line in f:
            values = line.split(",") #it's a CSV
            entry = dict(zip(headers, values))
            sample_name = entry["sample_Xtile_Ytile"]
            sample_id, xtile_id, ytile_id = sample_name.split("_")
            xtile_id = int(xtile_id[1:])
            ytile_id = int(ytile_id[1:])
            xx = int(entry["X.X"])
            yy = int(entry["Y.Y"])
            x = x_tile_size * (xtile_id-1) + xx
            y = y_tile_size * (ytile_id-1) + yy
            cell_type = cluster_annotations[entry["Imaging phenotype cluster ID"]]
            codex.add_cell(sample_id, x, y, cell_type)
    return codex

def read_colorectal_csv(colorectal_csv_path: str) -> CodexData:
    with open(colorectal_csv_path) as f:
        codex = CodexData()
        headers = f.readline().split(",")
        for line in f:
            values = line.split(",") #it's a CSV
            entry = dict(zip(headers, values))
            region = entry["Region"]
            x = int(entry["X:X"])
            y = int(entry["Y:Y"])
            cell_type = entry["ClusterName"]
            codex.add_cell(region, x, y, cell_type)
    return codex
