import pandas as pd
from re import split


def read_file(file, cog_pos):
    """
    read processed file and get only location of feature and assigned COG
    """
    data = pd.read_csv(file, sep='\t')
    attributes = [split(r"[=,;]", x) for x in list(data["attribute"])]
    new_data = data[["start", "end"]]
    new_data.insert(2, "COG", [og[cog_pos] for og in attributes], True)
    return [data, new_data]


def consensus(em_file, om_file, batch_file, fasta_file, get_pseudo=False, get_ncrna=False, gff_file=None, output_dir=""):
    """
    Improve the functional annotation of the bacterial genome using a consensus of three programs:
    eggNOG-mapper, Operon-mapper and Batch CD-Search. Function saves all predicted features and COG assignments
    and prepares the outputs file for visualization with DNAPlotter
    :param em_file: the path to Eggnog-mapper processed file
    :param om_file: the path to Operon-mapper processed file
    :param batch_file: the path to Batch CD-Search processed file
    :param fasta_file: the path to genomic sequence
    :param output_dir: the output file
    :type get_pseudo: bool
    :type get_ncrna: bool
    :param gff_file: the path to gff file where all features are stored
    :return:  file with functional annotation of the bacterial genome
    """
    # read eggnog_mapper file
    [em_data, new_em_data] = read_file(em_file, 3)
    # read operon_mapper file
    [om_data, new_om_data] = read_file(om_file, 3)
    # read batch cd search file
    [batch_data, new_batch_data] = read_file(batch_file, 1)

    # merge the three dataframes to save all predicted features into one dataframe
    # cog_x = eggnog-mapper #cog_y = operon-operon-mapper #cog = batch cd-search
    new_df = pd.merge(new_em_data, new_om_data, on="start", how="outer")
    new_df = (pd.merge(new_df, new_batch_data, on="start", how="outer")).fillna("-")

    # create new DataFrame to store the necessary information
    df = pd.DataFrame(columns=["seqname", "source", "type", "start", "end", "score", "strand", "frame", "attribute"])

    # iterate through features
    for row in new_df.index:
        # get assigned COGs, number of NaN and start of the feature
        cogs = [new_df["COG_x"][row], new_df["COG_y"][row], new_df["COG"][row]]
        nan = cogs.count("-")
        start = new_df["start"][row]
        # find out which tools match in the COG assignment
        idxs = [[cogs[:idx].index(item), idx] for idx, item in enumerate(cogs) if
                item in cogs[:idx]]

        # if all three tools have assigned the COG
        if nan == 0:
            # if all three tools match in assignment, add eggnog-mapper
            if len(idxs) != 1:
                df = pd.concat([df, em_data.loc[em_data.start == start, :]], ignore_index=True)

            # two tools match in assignment
            else:
                #  if eggnog and batch or operon -> add eggnog, if operon and batch -> add operon
                df = pd.concat([df, om_data.loc[om_data.start == start, :]], ignore_index=True) if idxs[0] == [1, 2] \
                    else pd.concat([df, em_data.loc[em_data.start == start, :]], ignore_index=True)

        # if only one tool has assigned the COG
        elif nan == 2:
            # find out which one was it and add it
            which = [number for number in [0, 1, 2] if number not in idxs[0]][0]
            programs = ['em_data', 'om_data', 'batch_data']
            df = pd.concat([df, (locals()[programs[which]]).loc[(locals()[programs[which]]).start == start, :]],
                           ignore_index=True)

        # if two tools have assigned the COG
        elif nan == 1:
            nan_position = cogs.index("-")
            # add eggnog-mapper or operon-mapper
            df = pd.concat([df, om_data.loc[om_data.start == start, :]], ignore_index=True) if nan_position == 0 else \
                pd.concat([df, em_data.loc[em_data.start == start, :]],
                          ignore_index=True)

        # no tool has assigned the COG, either it is not assigned or it is not a CDS -> add operon-mapper
        else:
            df = pd.concat([df, om_data.loc[om_data.start == start, :]], ignore_index=True)

    # add pseudogenes and/or ncRNA
    if get_pseudo or get_ncrna:
        df = get_features(gff_file, df, get_pseudo, get_ncrna)

    # save the created dataframe into new file and add genomic sequence
    df.to_csv(output_dir + '/file_to_plot.txt', sep='\t', index=False, header=False)
    fasta_data = (open(fasta_file, "r")).read()
    with open(output_dir + '/file_to_plot.txt', 'a') as my_file:
        my_file.write('\n' + fasta_data)
    my_file.close()


def get_features(gff_file, df, get_pseudo, get_ncrna):
    """
    change the feature type to a pseudogene according to information in gff_file
    and add ncRNA feature to the dataframe
    """
    gff_file = pd.read_csv(gff_file, sep="\t", header=None, comment="#", names=("seqname", "source", "type", "start",
                                                                                "end", "score", "strand", "frame",
                                                                                "attribute"))
    if get_pseudo:
        pseudogenes = gff_file.loc[gff_file['type'] == 'pseudogene']
        for row in pseudogenes.index:
            start = pseudogenes.start[row]
            df.loc[df['start'] == start, 'type'] = 'pseudogene'
    if get_ncrna:
        ncrnas = gff_file.loc[gff_file['type'] == 'ncRNA']
        df = pd.concat([df,ncrnas], ignore_index=True)

    return df
