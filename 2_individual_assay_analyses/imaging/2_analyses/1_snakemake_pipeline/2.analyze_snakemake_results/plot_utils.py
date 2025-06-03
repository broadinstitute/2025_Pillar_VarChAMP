from functools import reduce
import operator
import numpy as np
import polars as pl
import matplotlib.pyplot as plt
import seaborn as sns


## Letter dict to convert well position to img coordinates
letter_dict = {
    "A": "01",
    "B": "02",
    "C": "03",
    "D": "04",
    "E": "05",
    "F": "06",
    "G": "07",
    "H": "08",
    "I": "09",
    "J": "10",
    "K": "11",
    "L": "12",
    "M": "13",
    "N": "14",
    "O": "15",
    "P": "16",
}


## Channel dict to map channel to cellular compartments
channel_dict = {
    "DAPI": "1",
    "GFP": "2",
    "AGP": "3",
    "Mito": "4",
    "Brightfield1": "5",
    "Brightfield2": "6",
    "Brightfield": "7",
}


color_map = {
    'TC': 'slategrey', # Grey for controls
    'NC': 'gainsboro', 
    'PC': 'plum',
    'cPC': 'pink',
    'cNC': 'lightgrey',
    'allele': 'salmon',  # Tomato for disease
    'disease_wt': 'lightskyblue',  # Skyblue for reference
    '': 'white'  # White for missing wells
}


def plot_platemap(
    df,
    plate_name,
    well_pos_col="well_position",
    # this is the column to color by (categorical or continuous)
    value_col="node_type",
    # these columns will be concatenated into the annotation text
    label_cols=("gene_allele",),
    value_type="categorical",   # or "continuous"
    ax=None,
    continuous_cmap="vlag",  # matplotlib colormap for continuous mode
    categorical_colors=color_map,     # dict for categorical → color
    grid_square=None
):
    # 1) build empty 16×24 grid
    rows = list("ABCDEFGHIJKLMNOP")
    cols = [f"{i:02d}" for i in range(1,25)]
    plate_grid = (
        pl.DataFrame({c: [""]*16 for c in cols})
          .with_row_index("row_index")
          .unpivot(index="row_index", on=cols, variable_name="col_label", value_name="_")
          .with_columns(
              pl.col("row_index").map_elements(lambda i: rows[i], return_dtype=pl.Utf8).alias("row_label")
          )
    )
    # display(plate_grid)

    # 2) extract row/col from your df’s well_position
    df2 = df.with_columns([
        pl.col(well_pos_col).str.head(1).alias("row_label"),
        pl.col(well_pos_col).str.slice(1).alias("col_label")
    ])

    # 3) join
    plate = plate_grid.join(df2, on=["row_label","col_label"], how="left")

    # 4) pivot out two matrices:
    #    A) data matrix for coloring
    #    B) text matrix for annotation
    # first build annotation text by concatenating label_cols
    plate = plate.with_columns(
        reduce(
            lambda acc, c: acc + "\n" + \
            pl.col(c).round(2).cast(pl.Utf8).fill_null(""),
            label_cols[1:],
            pl.col(label_cols[0]).fill_null("").str.replace_all("_", "\n")
        ).alias("_annot")
    )
    # display(plate)

    # pivot color‐matrix
    data_matrix = plate.pivot(
        index="row_label", on="col_label", values=value_col
    )

    # pivot annotation‐matrix
    annot_matrix = plate.pivot(
        index="row_label", on="col_label", values="_annot"
    ).fill_null("")

    # convert to numpy
    # drop the implicit “row_label” column in position 0
    data = data_matrix[:,1:].to_numpy()
    ann = annot_matrix[:,1:].to_numpy()

    # 5) choose coloring
    if value_type == "categorical":
        if categorical_colors is None:
            raise ValueError("Must supply categorical_colors when value_type='categorical'")
        # map each category in data to its color
        # build vectorized map
        cmap_array = np.vectorize(lambda x: categorical_colors.get(x, "white"))(data)
        # For seaborn we draw a dummy zero‐matrix
        plot_data = np.zeros_like(data, dtype=float)
        cmap = None
    else:
        # continuous: data is numeric
        plot_data = data.astype(float)
        cmap = continuous_cmap
        cmap_array = None

    # 6) plot
    if ax is None:
        fig, ax = plt.subplots(1,1,figsize=(35,14))
        
    sns.heatmap(
        plot_data,
        ax=ax,
        annot=ann if value_type=="categorical" else None,
        fmt="",
        cmap=cmap,
        cbar=(value_type=="continuous"),
        # linewidths=0,
        # linecolor="white",
        square=True,
        annot_kws={"size":9, "color": "black"}
    )

    # if categorical: overlay colored rectangles
    if value_type=="categorical":
        for i in range(cmap_array.shape[0]):
            for j in range(cmap_array.shape[1]):
                ax.add_patch(plt.Rectangle(
                    (j, i), 1, 1,
                    color=cmap_array[i,j],
                    # ec="black"
                ))
    else:
        # create combined annotation: value + other labels
        # you could easily extend to show gene_allele too by rebuilding ann
        for i in range(ann.shape[0]):
            for j in range(ann.shape[1]):
                txt = ann[i,j]
                # if you want gene_allele too: append "\n"+ann[i,j]
                ax.text(
                    j+0.5, i+0.5, txt,
                    ha="center", va="center", fontsize=9.5, color="black"
                )

    if grid_square is not None:
        grid_sq_mat = plate.pivot(
            index="row_label", on="col_label", values=grid_square
        )[:,1:]#.to_numpy()
        for i in range(grid_sq_mat.shape[0]):
            for j in range(grid_sq_mat.shape[1]):
                if grid_sq_mat[i,j] is not None and grid_sq_mat[i,j]>=1:
                    ax.add_patch(plt.Rectangle(
                        (j, i), 1, 1,
                        linewidth=2, edgecolor="red", facecolor="none"
                    ))

    # 7) finalize axes
    ax.set_title(f"384-Well Plate: {plate_name}")
    ax.set_xlabel("Column")
    ax.set_ylabel("Row")
    ax.set_xticks(np.arange(len(cols))+0.5)
    ax.set_xticklabels(cols, rotation=0)
    ax.set_yticks(np.arange(len(rows))+0.5)
    ax.set_yticklabels(rows, rotation=0)
    # plt.tight_layout()
    # plt.show()
    return plate