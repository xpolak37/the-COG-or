import sys
import getopt
import program_processor
import track_manager
import consensus
import os


def cogor():
    organism_name = None
    input_dir = os.getcwd()
    output_dir = os.getcwd()
    manager = False

    argv = sys.argv[1:]

    try:
        opts, args = getopt.getopt(argv,"n:i:o:t")

    except:
        print("Something wrong with the arguments")
        sys.exit(2)

    try:
        for opt, arg in opts:
            if opt in ["-n"]:
                organism_name = arg
            elif opt in ["-i"]:
                input_dir = arg
            elif opt in ["-o"]:
                output_dir = arg
            elif opt in ["-t"]:
                if arg in ["true", "yes", "True", "t", ""]:
                    manager = True

    except:
        print("Something wrong with the arguments")
        sys.exit(2)

    try:
        # Program processor
        program_processor.em_processor(organism_name,input_dir + "/" + organism_name + "_eggnog.gff", input_dir + "/" +
                                       organism_name + "_cds.txt", output_dir)
        program_processor.om_processor(organism_name, input_dir + "/" + organism_name + "_orf_operon.txt", input_dir +
                                       "/" + organism_name + "_cog_operon.txt", output_dir)
        program_processor.batch_processor(organism_name, input_dir + "/" + organism_name + "_batch.txt", output_dir)

        # Consensus
        consensus.consensus(output_dir + "/em_" + organism_name + ".gff",
                            output_dir + "/om_" + organism_name + ".gff",
                            output_dir + "/batch_" + organism_name + ".gff",
                            input_dir + "/" + organism_name + ".fasta",
                            get_pseudo=True, get_ncrna=True,
                            gff_file=input_dir + "/" + organism_name + ".gff3",
                            output_dir=output_dir)

        # Track manager
        if manager:
            track_manager.get_track_template(output_dir=output_dir)
            track_manager.get_legend(output_dir=output_dir)

    except:
        print("Something wrong with your files.")
        sys.exit(2)


cogor()
