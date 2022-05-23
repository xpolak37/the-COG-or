import seaborn as sns
import pkg_resources
from PIL import Image, ImageDraw, ImageFont
import os


def get_track_template(pos_track=(0.95, 0.90, 0.85, 0.80), size=10.0, output_dir=os.getcwd()):
    """
    Generate file for Track Manager option in DNAPlotter
    :param pos_track: the positions for plotting features (CDS forward strand, CDS reverse strand, pseudogenes, RNA genes)
    :param size: the size of track
    :param output_dir: the output file
    :return: track template file
    """
    # position vector: 27x for CDS forward strand, 27x for CDS reverse strand, 27x for pseudogenes, 3x for RNA genes
    pos = [str(pos_track[0]), str(pos_track[1])] * 27
    pos.extend([str(pos_track[2])] * 27)
    pos.extend([str(pos_track[3])] * 3)

    # size vector
    size = [str(size)] * 84

    # forward strand: CDS...27x true, 27x false, pseudogenes...27x true, RNAs...3x true
    fwd = ["true", "false"] * 27
    fwd.extend(["true"] * 27)
    fwd.extend(["true"] * 3)

    # reverse strand: CDS...27x false, 27x true, pseudogenes...27x true, RNAs...3x true
    rev = ["false", "true"] * 27
    rev.extend(['true'] * 27)
    rev.extend(["true"] * 3)

    NOT, ANY = ["false"] * 84, ["false"] * 84

    # features to plot: CDS, pseudogene, tRNA, rRNA, ncRNA
    KEY = ["CDS"] * 54
    KEY.extend(['pseudogene'] * 27)
    KEY.extend(["tRNA", "rRNA", "ncRNA"])

    # qualifier to distinguish features: category for CDS and pseudogenes, null for RNAs
    QUAL = ["CAT"] * 81
    QUAL.extend(["null"] * 3)

    # 26 COG categories, '-' for COG unknown, 'null' for RNAs
    VAL = ["J", "J", "A", "A", "K", "K", "L", "L", "B", "B", "D", "D", "Y", "Y", "V", "V", "T", "T", "M", "M", "N", "N",
           "Z", "Z","W", "W", "U", "U", "O", "O", "X", "X", "C", "C", "G", "G", "E", "E", "F", "F", "H", "H", "I", "I",
           "P", "P", "Q", "Q", "R", "R", "S", "S", "-", "-",
           "J", "A", "K", "L", "B", "D", "Y", "V", "T", "M", "N", "Z", "W", "U", "O", "X", "C", "G", "E", "F", "H", "I",
           "P", "Q", "R", "S", "-"]
    VAL.extend(["null"] * 3)

    # color palette: 20 colors have been chosen from seaborn palette and other ten were selected subjectively
    palette = sns.color_palette('tab20')
    palette2 = [(255, 209, 49), (244, 172, 50), (219, 179, 177), (165, 190, 0), (61, 82, 213), (255, 238, 221),
                (217, 3, 104), (250, 243, 62), (30, 252, 30)]
    palette = [str(round(color[0] * 250)) + ":" + str(round(color[1] * 250)) + ":" + str(round(color[2] * 250))
               for color in palette]
    palette2 = [str(color[0]) + ":" + str(color[1]) + ":" + str(color[2]) for color in palette2]
    palette = palette + palette2

    # add color palette to color vector
    COL = [0] * 78
    for i in range(2):
        COL[i:52:2] = palette[:26]
    COL[53:80] = palette[:26]
    COL[52:53] = ["0:0:0", "0:0:0"]
    COL.extend(["0:0:0", palette[26], palette[27], palette[28]])

    # join all vectors into string
    string = ''
    for i in range(84):
        string = string + pos[i] + "\t" + size[i] + "\t" + fwd[i] + "\t" + rev[i] + "\t" \
                 + NOT[i] + "\t" + ANY[i] + "\t" + KEY[i] + "\t" + QUAL[i] + "\t" + VAL[i] + \
                 "\t" + str(COL[i]) + "\n"

    with open(output_dir + "/track_template", "w") as file_to_save:
        file_to_save.write(string)
        file_to_save.close()


def get_legend(output_dir=os.getcwd()):
    """
    Create a legend for the genome map
    """
    # color palette: 20 colors have been chosen from seaborn palette and other ten were selected subjectively
    palette = sns.color_palette('tab20')
    palette2 = [(255, 209, 49), (244, 172, 50), (219, 179, 177), (165, 190, 0), (61, 82, 213), (255, 238, 221),
                (0, 0, 0), (217, 3, 104), (250, 243, 62), (30, 252, 30)]
    palette = [((round(tuple[0] * 250)), round(tuple[1] * 250), round(tuple[2] * 250)) for tuple in palette]
    palette = palette + palette2

    # import description for the legend
    file = pkg_resources.resource_filename(__name__, 'COGor-data/legend_text.csv')
    CATS = (open(file, "r")).read()
    CATS = CATS.split('\n')
    # create a white image
    img = Image.new("RGB", (1700, 2500), "white")
    image_edit = ImageDraw.Draw(img)
    myFont = ImageFont.truetype('arial.ttf', 50)

    # add individual objects to legend
    start = 50
    for ind in range(len(palette)):
        image_edit.rectangle((50, start, 50 + 80, start + 80), fill=(palette[ind][0], palette[ind][1], palette[ind][2]))
        image_edit.text((150, start + 15), CATS[ind], font=myFont, fill=(0, 0, 0))
        start = start + 80
    img.save(output_dir + '/legend.jpg')
