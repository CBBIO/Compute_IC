# ----------------------------------------------------------------------------
# Textual and Graphical Information Content from GO_IDs
# Created By  : Israel Barrios, Ana Rojas Mendoza
# Credit: CBBio Lab (Centro Andaluz de Biología del Desarrollo)
# Contact: israel.barrios@csic.es
# Created Date: 4/5/2023
# Status: Prototype
# version ='0.2'
# ---------------------------------------------------------------------------
# Klopfenstein DV, Zhang L, Pedersen BS, ... Tang H GOATOOLS: A Python library for Gene Ontology analyses Scientific reports | (2018) 8:10872 | DOI:10.1038/s41598-018-28948-z
# Martin Larralde; Philipp A.; Alex Henrie; Daniel Himmelstein; Dave Lawrence; Rafal Wojdyla; Spencer Mitchell; Tatsuya Sakaguchi althonos/pronto: v2.5.4 | 10.5281/zenodo.7814219
# Pauli Virtanen, Ralf Gommers, Travis E. Oliphant, Matt Haberland, Tyler Reddy, David Cournapeau, Evgeni Burovski, Pearu Peterson, Warren Weckesser, Jonathan Bright, Stéfan J. van der Walt, Matthew Brett, Joshua Wilson, K. Jarrod Millman, Nikolay Mayorov, Andrew R. J. Nelson, Eric Jones, Robert Kern, Eric Larson, CJ Carey, İlhan Polat, Yu Feng, Eric W. Moore, Jake VanderPlas, Denis Laxalde, Josef Perktold, Robert Cimrman, Ian Henriksen, E.A. Quintero, Charles R Harris, Anne M. Archibald, Antônio H. Ribeiro, Fabian Pedregosa, Paul van Mulbregt, and SciPy 1.0 Contributors. (2020) SciPy 1.0: Fundamental Algorithms for Scientific Computing in Python. Nature Methods, 17(3), 261-272.

import os
import argparse
import IC_lib as IC
if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser(
        prog="GOIC", description="GO Information Content Calculator"
    )
    parser.add_argument(
        "inputfile", help="Tab separated file with \"ID\tGO_ID\" rows")
    parser.add_argument("-o","--outputpath",help="Output file.\n A pickled dict if --precompute is selected or a tsv with ICs if it's not",default=None)
    parser.add_argument("--precompute",help="Precompute IC values into PATH",action="store_true")
    parser.add_argument("--precomputed_db",help="Path to pickled dict with precomputed IC values",default=None)
    args = parser.parse_args()

    if (args.precompute):
        if (args.outputpath is not None) and (os.path.exists(args.inputfile)):
            outputpath=args.outputpath
            if os.path.isdir(outputpath):
                os.mkdir(outputpath)
                outputpath=outputpath+"/precomputed_IC_file.pickle"
            IC.precompute_data(args.inputfile,outputpath)
        else:
            print("""
            An universe of GOs need to be used to precompute ICs.\n
            Use 'python compute_IC.py --precompute ANNOTATIONSFILEPATH OUTPUTPATH'\n
            ANNOTATIONSFILEPATH -> Can be a \\n separated list of GOs or a GAF file
            """)
            exit()
    else:
        if args.precomputed_db is not None:    
            if os.path.exists(args.inputfile) and os.path.exists(args.precomputed_db):
                print("Using custom precomputed DB")
                IC.precalc_IC(args.inputfile,args.precomputed_db)
                
            else:
                if not os.path.exists(args.inputfile):
                    print("Wrong inputfile")
                elif not os.path.exists(args.precomputed_db):
                    print("""
                    Wrong precomputed pickle\n
                    If you want to create it, use 'python compute_IC.py --precompute ANNOTATIONSFILEPATH OUTPUTPATH'
                    """)
        else:
            print("Using default GOA_UNIPROT precomputed DB")
            IC.precalc_IC(args.inputfile)