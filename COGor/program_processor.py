from re import search, split
import pandas as pd
from Bio import SeqIO
import os
import pkg_resources
import re


def em_processor(organism_name, em_file, cds_file, output_dir=os.getcwd()):
    """
    Process the output file (decorated.gff) from eggNOG-mapper tool into more structured COGor-data.
    The outputs of this function is file in gff format that contains a suitable header with information about CDSs with
    assigned COG by eggNOG-mapper
    :type organism_name: str
    :param em_file: the path to eggNOG-mapper output file
    :param cds_file: the path to eggNOG-mapper input file
    :param output_dir: the output directory
    :return: processed file
    """
    em_data = pd.read_csv(em_file, sep="\t", header=None, comment="#", names=("seqname", "source", "type", "start",
                                                                              "end", "score", "strand", "frame",
                                                                              "attribute"))
    cds_data = SeqIO.parse(cds_file, "fasta")
    # get only header of each CDS
    desc = [record.description for record in cds_data]
    # iterate through CDSs header
    cogs_file = pkg_resources.resource_filename(__name__, 'COGor-data/cogs.txt')
    cat_data = (open(cogs_file, "r")).readlines()

    for record in desc:
        seq_id = record[0:record.index(" [")]
        # include all possible locations forward/reverse strand and join
        dic = {'pattern': "complement\((.*?)\)", 'strand': "-"} if "complement" in record else \
            {'pattern': "location=(.*?)]", 'strand': "+"}
        location = search(dic['pattern'], record).group(1)

        if "<" in location:
            location = location[1:]

        if "join" in location:
            location = location[:-1] if location[-1] == ")" else location
            location = (location[5:].split(".."))
            location = [location[0], location[len(location) - 1]]
        else:
            location = location.split("..")

        if ">" in location[1]:
            location[1] = location[1][1:]

        # add the information about location to the corresponding row of eggnog-mappers outputs
        em_data.loc[em_data.seqname == seq_id, ["seqname", "strand", "start", "end"]] = \
            [seq_id[seq_id.index("|") + 1:seq_id.index("_prot")], dic['strand'], location[0], location[1]]

    for row in em_data.index:
        # get only useful information about each CDS: feature_id, name, COG, COG category
        attribute = split(r"[=,;]", em_data["attribute"][row])
        seq_id = attribute[attribute.index("ID") + 1]
        dic = {'COG': (attribute[attribute.index("em_OGs") + 1].split("@"))[0],
               'desc': attribute[attribute.index("em_desc") + 1],
               'cat': attribute[attribute.index("em_COG_cat") + 1]}
        # if cog is from COG database
        try:
            index = [i for i, s in enumerate(cat_data) if dic["COG"] in s][0]
            end_index = re.search('\n', cat_data[index]).start()
            #dic["cat"] = cat_data[index][8:end_index]
            dic["cat"] = cat_data[index][8]
        # else - it is from eggNOG or ROG
        except IndexError:
            dic["cat"] = dic["cat"][0]
        try:
            dic['name'] = attribute[attribute.index("em_Preferred_name") + 1]
        except ValueError:
            dic['name'] = '-'

        dic['cat'] = 'S' if dic['cat'] == 'None' else dic['cat']
        em_data.loc[row, "attribute"] = "".join(["ID=", seq_id, ";COG=", dic['COG'], ";CAT=", dic['cat'],
                                                 ";name=", dic['name'], ";desc=", dic['desc']])

    return em_data.to_csv(output_dir + '/em_' + organism_name + '.gff', sep='\t', index=False)


def om_processor(organism_name, orf_file, cog_file, output_dir=os.getcwd()):
    """
    Process the outputs files (ORF_coordinates.txt and predicted_COGs.txt) from Operon-mapper into more structured COGor-data.
    The outputs of this function is file in gff format that contains a suitable header with information about all
    predicted features
    :type organism_name: str
    :param orf_file: the path to Operon-mapper outputs file ORFs_coordinates.txt
    :param cog_file: the path to Operon-mapper outputs file predicted_COGs.txt
    :param output_dir: the output directory
    :return: processed file
    """
    orf_data = pd.read_csv(orf_file, sep="\t", header=None, comment="#", names=("seqname", "source", "type", "start",
                                                                                "end", "score", "strand", "frame",
                                                                                "attribute"))
    cog_data = pd.read_csv(cog_file, sep="\t", header=None, comment="#", names=("ID", "COG", "category"))
    cogs_file = pkg_resources.resource_filename(__name__, 'COGor-data/cogs.txt')
    cat_data = (open(cogs_file, "r")).readlines()

    for row in orf_data.index:
        # iterate through all features in ORF file and save the relevant information from the COG file
        attribute = split(r"[=,;]", orf_data["attribute"][row])
        feature_id = attribute[attribute.index("ID") + 1]
        try:
            cog = cog_data.loc[cog_data["ID"] == feature_id, "COG"].values[0]

            # if cog is from COG database
            try:
                index = [i for i, s in enumerate(cat_data) if cog in s][0]
                end_index = re.search('\n', cat_data[index]).start()
                #CAT = cat_data[index][8:end_index]
                CAT = cat_data[index][8]

            # else - it is from eggNOG or ROG
            except IndexError:
                CAT = (cog_data.loc[cog_data["ID"] == feature_id, "category"].values[0])[1][0]

            dic = {'cat': 'S', 'desc': '-'} if "ROG" in cog else \
                {'cat': CAT,'desc': cog_data.loc[cog_data["ID"] == feature_id, "category"].values[0].split("] ")[1]}

            orf_data.loc[row, "attribute"] = "".join(
                ["ID=", feature_id, ";COG=", cog, ";CAT=", dic['cat'], ";desc=", dic['desc']])

        except:
            orf_data.loc[row, "attribute"] = "".join(["ID=", feature_id, ";COG=", "-", ";CAT=", "-"])

    return orf_data.to_csv(output_dir + '/om_' + organism_name + '.gff', sep='\t', index=False)


def batch_splitter(organism_name, gene_file,output_dir=os.getcwd()):
    """
    :type organism_name: str
    :param gene_file: the path to file to be split
    :param output_dir: the output directory
    :return: two split files
    """
    records = list(SeqIO.parse(gene_file, "fasta"))
    if len(records) > 4000:
        SeqIO.write(records[0:3500], organism_name + "_genes1.fasta", "fasta")
        SeqIO.write(records[3500:len(records)], output_dir + organism_name + "_genes2.fasta", "fasta")
    else:
        print("The file does not need to be split because it does not contain more than 4000 sequences.")


def batch_merger(organism_name, file1, file2,output_dir=os.getcwd()):
    """
    :type organism_name: str
    :param file1: the first annotated file from Batch-CD Search
    :param file2: the second annotated file from Batch-CD Search
    :param output_dir: the output directory
    :return: a merged file
    """
    data1, data2 = open(file1).read(), open(file2).read()
    data1 = data1[data1.index("Q#1"):len(data1)]
    data2 = data2[data2.index("Q#1"):len(data2)]
    data1 += data2

    with open(output_dir + organism_name + "_merged_hitdata.txt", "w") as file:
        file.write(data1)
        file.close()


def batch_processor(organism_name, batch_file, output_dir=os.getcwd()):
    """
    Process the outputs file (hitdata.txt) from Batch CD-Search tool into more structured COGor-data.
    The outputs of this function is file in gff format that contains a suitable header with information about CDSs with
    assigned COG by Batch CD-Search
    :type organism_name: str
    :param batch_file: the path to Batch CD-Search outputs file hitdata.txt
    :param output_dir: the output file
    :return: processed file
    """

    # create new DataFrame to store the necessary information
    batch_gff = pd.DataFrame(
        columns=["seqname", "source", "type", "start", "end", "score", "strand", "frame", "attribute"])
    batch_data = (open(batch_file).read())
    batch_data = (batch_data[batch_data.index("Q#"):len(batch_data) - 1]).split('\n')
    query = ''
    cogs_file = pkg_resources.resource_filename(__name__, 'COGor-data/cogs.txt')
    cogs_data = (open(cogs_file, "r")).readlines()

    # iterate through queries
    for row in batch_data:
        new_query = search('Q#\d+', row).group(0)

        # if COG is not assigned or new_line is duplicate, continue to next iteration
        if not ('\tspecific\t' in row) or query == new_query:
            continue

        query = new_query
        seq_id = row[row.index("|") + 1:row.index("_prot")]
        id = row[row.index("|") + 1:row.index(" [")]
        id = "".join(["ID=", id])

        # include all possible locations forward/reverse strand and join
        dic = {'pattern': "complement\((.*?)\)", 'strand': '-'} if 'complement' in row else \
            {'pattern': "location=(.*?)]", 'strand': '+'}
        location = search(dic['pattern'], row).group(1)

        if "<" in location:
            location = location[1:]

        if "join" in location:
            location = location[:-1] if location[-1] == ")" else location
            location = (location[5:].split(".."))
            location = [location[0], location[len(location) - 1]]

        else:
            location = location.split("..")
            location = [location[0], location[1]]

        if ">" in location[1]:
            location[1] = location[1][1:]

        try:
            COG = "".join(["COG=", search("(COG\d+)", row).group(1)])

        except AttributeError:
            COG = "".join(["COG=", "-"])

        # update in COG 2021
        if COG == "COG=COG3512":
            COG = "COG=COG1343"

        index = [i for i, s in enumerate(cogs_data) if COG[4:] in s][0]
        end_index = re.search('\n', cogs_data[index]).start()
        #CAT = cogs_data[index][8:end_index]
        CAT = cogs_data[index][8]
        CAT = "".join(["CAT=", CAT])

        attribute = id + ";" + COG + ";" + CAT
        new_row = pd.DataFrame({"seqname": [seq_id], "source": ["unknown"], "type": ["CDS"], "start": [location[0]], "end": [location[1]],
                   "score": ["."], "strand": [dic["strand"]], "frame": ["0"], "attribute": [attribute]})

        batch_gff = pd.concat([batch_gff,new_row], ignore_index=True)
    return batch_gff.to_csv(output_dir + '/batch_' + organism_name + '.gff', sep='\t', index=False)
