#EXAMPLE

from matplotlib import pyplot as plt
import read_csvs

#data extracted from Table S2 and put into a CSV file for easier use
SPLEEN_ANNOTATIONS_PATH = r"data\spleen\cluster_annotations.csv"
#raw data straight from the paper
RAW_SPLEEN_CSV_PATH = r"data\spleen\raw_codex_data.csv"

codex = read_csvs.read_raw_spleen_csv(RAW_SPLEEN_CSV_PATH, SPLEEN_ANNOTATIONS_PATH)

sample = codex.samples["BALBc-3"].get_subsample(2000, 2000, 4000, 4000)

for phenotype in sample.get_top_phenotypes(): #get_top_phenotypes so the legend is sorted by frequency
    xs, ys = sample.get_cells(phenotype)
    plt.plot(xs, ys, ".", label=phenotype)

plt.gca().set_aspect("equal")
plt.legend(bbox_to_anchor=(1.04, 1)) #put the legend outside the plot
plt.show()