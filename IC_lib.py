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
import re
import pickle
import gzip
from functools import lru_cache
import math


from pronto import Ontology
import matplotlib.pyplot as plt
import requests
import numpy
from scipy import stats
import p_tqdm,tqdm

# Not requiring gui so we select agg backend
import matplotlib
matplotlib.use("Agg")

###############################################################################
# CONFIG CONSTANTS
###############################################################################
BASIC_OBO = False  # True for using go-basic.obo, False for go.obo
CPU_COUNT = int(os.cpu_count())  # Change this to select the number of cores


###############################################################################
# UTILS
###############################################################################
# Creates a new folder if it doesn't exists
def new_folder(FOLDER):
    if not os.path.exists(FOLDER):
        os.mkdir(FOLDER)
    return FOLDER

# Remote file requests


def get_remote_file(URL, FILEPATH):
    if not os.path.exists(FILEPATH):
        data = requests.get(URL).text
        with open(FILEPATH, "w") as fwrite:
            fwrite.write(data)
    return FILEPATH

# Pickle object to file


def pickle_object(data_object, FILEPATH):
    with open(FILEPATH, "wb") as fwrite:
        pickle.dump(data_object, fwrite)

# Unpickle from file


def unpickle_object(FILEPATH):
    with open(FILEPATH, "rb") as fread:
        data = pickle.load(fread)
    return data


###############################################################################
# OBO
###############################################################################
# Simple OBO parser, Only GO_ids are returned
def parse_obo(FILEPATH):
    regex = re.compile("alt_id: (GO:.*)\n")
    with open(FILEPATH, "r") as fread:
        data = [
            [y]
            + [
                y.split(": ")[1].strip()
                for y in x.strip().split("\n")
                if y and ":" in y
            ][1:3]
            for x in fread.read().replace(";", ",").split("\n\n")
            if x and ("[Term]" in x)
            for y in (
                [x.strip().split("\n")[1].split(": ")[
                    1].strip()] + regex.findall(x)
            )
        ]
    return {x[0]: x[1::] for x in data}


# Get term depth
# A lite convoluted, it needs work
@lru_cache(maxsize=None)
def get_depth(pronto_Term,aspect):
    entity = pronto_Term
    depth = []
    i = 0
    entity_iterator = [entity]

    while entity_iterator:
        entity_aux = set()
        for x in entity_iterator:
            parents = [y for y in x.superclasses(distance=1, with_self=False) if annotations[aspect].get(y.id,False)]
            entity_aux.update(parents)
        if entity_aux:
            depth.append(entity_aux)
        entity_iterator = entity_aux
        i = i + 1

    return len(depth)

###############################################################################
# ANNOTATIONS
###############################################################################

# The output of this function IS NOT A SET so IT HAS DUPLICATES


def read_gpad(FILEPATH):
    with gzip.open(FILEPATH, "rt") as fread:
        try:
            data = fread.read()
        except OSError:
            print("Annot file is not a valid gzip file by OSError")
            print("Trying to read as a text file")
            with open(FILEPATH, "r") as fread_txt:
                data = fread_txt.read()
    gpad_go = [
        x.split("\t")[3] for x in data.split("\n") if x and (not x.startswith("!"))
    ]
    return gpad_go


def compile_universe(all_annot):
    universe = {
        "biological_process": {},
        "cellular_component": {},
        "molecular_function": {},
    }
    for goterm in all_annot:
        if obo.get(goterm, None) is not None:
            aspect = obo[goterm][1]
            universe[aspect][goterm] = universe[aspect].get(goterm, 0) + 1
    print("Compiling annotations")
    print("Total of annotations readed", len(all_annot))
    for x in universe:
        universe[x]["total"] = sum(
            [universe[x][y] for y in universe[x] if y != "total"]
        )
        print("", x, universe[x]["total"])

    return universe


def get_universe(progress=False):
    ANNOTATIONS = ANNOTATIONS_FOLDER
    universe_PICKLE = UNIVERSE_COUNT_PICKLE
    abc2gpad = {
        "hsa": "goa_human.gpad.gz",  # human
        "mmu": "mgi.gpad.gz",  # mouse
        "dme": "fb.gpad.gz",  # fly
        "cgd": "cgd.gpad.gz",
        "dictybase": "dictybase.gpad.gz",
        "ecocyc": "ecocyc.gpad.gz",
        "fb": "fb.gpad.gz",
        "genedb_lmajor": "genedb_lmajor.gpad.gz",
        "genedb_pfalciparum": "genedb_pfalciparum.gpad.gz",
        "genedb_tbrucei": "genedb_tbrucei.gpad.gz",
        "goa_chicken": "goa_chicken.gpad.gz",
        "goa_chicken_complex": "goa_chicken_complex.gpad.gz",
        "goa_chicken_isoform": "goa_chicken_isoform.gpad.gz",
        "goa_chicken_rna": "goa_chicken_rna.gpad.gz",
        "goa_cow": "goa_cow.gpad.gz",
        "goa_cow_complex": "goa_cow_complex.gpad.gz",
        "goa_cow_isoform": "goa_cow_isoform.gpad.gz",
        "goa_cow_rna": "goa_cow_rna.gpad.gz",
        "goa_dog": "goa_dog.gpad.gz",
        "goa_dog_complex": "goa_dog_complex.gpad.gz",
        "goa_dog_isoform": "goa_dog_isoform.gpad.gz",
        "goa_dog_rna": "goa_dog_rna.gpad.gz",
        "goa_human": "goa_human.gpad.gz",
        "goa_human_complex": "goa_human_complex.gpad.gz",
        "goa_human_isoform": "goa_human_isoform.gpad.gz",
        "goa_human_rna": "goa_human_rna.gpad.gz",
        "goa_pig": "goa_pig.gpad.gz",
        "goa_pig_complex": "goa_pig_complex.gpad.gz",
        "goa_pig_isoform": "goa_pig_isoform.gpad.gz",
        "goa_pig_rna": "goa_pig_rna.gpad.gz",
        "goa_uniprot_all_noiea": "goa_uniprot_all_noiea.gpad.gz",
        "japonicusdb": "japonicusdb.gpad.gz",
        "mgi": "mgi.gpad.gz",
        "pombase": "pombase.gpad.gz",
        "pseudocap": "pseudocap.gpad.gz",
        "reactome": "reactome.gpad.gz",
        "rgd": "rgd.gpad.gz",
        "sgd": "sgd.gpad.gz",
        "sgn": "sgn.gpad.gz",
        "tair": "tair.gpad.gz",
        "wb": "wb.gpad.gz",
        "xenbase": "xenbase.gpad.gz",
        "zfin": "zfin.gpad.gz",
    }
    if not os.path.exists(UNIVERSE_COUNT_PICKLE):
        print("DOWNLOADING ANNOTATION FILEs")
        all_annot = []
        for i, file in enumerate(abc2gpad):
            filepath = ANNOTATIONS + abc2gpad[file]
            if progress:
                print(filepath, str(i) + "/" + str(len(abc2gpad)))
            if not os.path.exists(filepath):
                wget_gz = "http://current.geneontology.org/annotations/{GZ}".format(
                    GZ=abc2gpad[file]
                )
                r = requests.get(wget_gz, stream=True)
                with open(filepath, "wb") as f:
                    for chunk in r.raw.stream(1024, decode_content=False):
                        if chunk:
                            f.write(chunk)
            gpad_go = read_gpad(filepath)
            all_annot.extend(gpad_go)
        universe = compile_universe(all_annot)
        with open(universe_PICKLE, "wb") as fwrite:
            pickle.dump(universe, fwrite)
    else:
        with open(universe_PICKLE, "rb") as fread:
            universe = pickle.load(fread)
    return universe


###############################################################################
# GENERAL
###############################################################################


def process_file(file_data):
    data = [x + [obo[x[1]][1], obo[x[1]][0]] for x in file_data if obo.get(x[1],None)]
    return data

def read_input(input_filepath):
    with open(input_filepath, "r") as fread:
        input_data = process_file(
            [x.split("\t")[0:2] for x in fread.read().split("\n") if x and not x.startswith("#")] 
        )
    return input_data

###############################################################################
# IC CALCULATIONs
###############################################################################
@lru_cache(maxsize=None)
def calculate_IC(x, aspect):
    universe = annotations[aspect]

    if godag.get(x):
        alt_ids = godag[x].alternate_ids | set([x])
        alt_ics = []
        for x in alt_ids:
            if universe.get(x, None) is not None:
                offprings = godag[x].subclasses(with_self=False)
                result = -1 * math.log2(
                    (universe[x] + sum([universe.get(y.id, 0)
                     for y in offprings]))
                    / (universe["total"])
                )
                alt_ics.append(result)
        if alt_ics:
            if min(alt_ics):
                return min(alt_ics)
            else:
                return 0.0
        else:
            return None
    else:
        return None


def crow_compute(row_data):
    prot_id, goterm, aspect, desc = row_data
    if godag.get(goterm, None):
        entity = godag[goterm]
        ic_row = [
            prot_id,
            aspect,
            goterm,
            len(sorted([sc for sc in entity.subclasses(with_self=False) if annotations[aspect].get(sc.id,False)])),
            get_depth(entity,aspect),
            calculate_IC(goterm, aspect),
            desc,
        ]
        return ic_row
    else:
        return None

def crow_precalc(row_data):
    prot_id, goterm, aspect, desc = row_data
    if godag.get(goterm, None):
        entity = godag[goterm]
        ic_row = [
            prot_id,
            aspect,
            goterm,
            precomputed_ics.get(goterm,["None"]*3)[0],
            precomputed_ics.get(goterm,["None"]*3)[1],
            precomputed_ics.get(goterm,["None"]*3)[2],
            desc,
        ]
        return ic_row
    else:
        return None



def compute_compute(input_data, annotations, godag):
    input_data_ic = [x for x in p_tqdm.p_umap(
        crow_compute, input_data, num_cpus=CPU_COUNT) if x]
    # input_data_ic = [crow_compute(x) for x in tqdm.tqdm(
    #     input_data) if x]
    return input_data_ic



def compute_precalc(input_data, annotations, godag):
    input_data_ic = [crow_precalc(x) for x in tqdm.tqdm(input_data) if x]
    return input_data_ic


def dump_ic_data(ic_data, outputfile="ic_data.tsv"):
    # Same format as GOATOOLS
    with open(TSVOUTPUTS_FOLDER+outputfile, "w") as fwrite:
        fwrite.write("\n".join(["\t".join([str(y) for y in x])
                     for x in ic_data if x])+"\n")

def plot_density(ic_data, filename):
    bins = numpy.linspace(0, 25, 1000)
    fig, ax = plt.subplots()
    ic_vector = {}
    go_set = set()
    aspect_color = {
        "biological_process": "red",
        "cellular_component": "blue",
        "molecular_function": "green"
    }
    for x in ic_data:
        if x and x[5] != "None" and x[2] not in go_set:
            ic_vector[x[1]] = ic_vector.get(x[1], []) + [float(x[5])]
            go_set.update([x[2]])

    plot_vector = []

    for aspect in ic_vector:
        if len(ic_vector[aspect]) > 1:
            n, x, _ = plt.hist(
                ic_vector[aspect],
                bins=bins,
                density=True,
                histtype="step",
                color=aspect_color[aspect],
                label=aspect,
            )
            density = stats.gaussian_kde(ic_vector[aspect])
            plot_vector.append([x, density, aspect])
            plt.clf()
    plt.title("Information Content " + filename)
    for x in plot_vector:
        plt.plot(
            x[0],
            x[1](x[0]),
            color=aspect_color[x[2]],
            label=x[2],
        )
    fig.legend()
    plt.xlabel("IC ( $-log_{2}(p(t))$ )")
    plt.savefig(PLOTS_FOLDER + "InformationContent_" + filename + ".jpeg")
    plt.close()


def precalc_IC(filepath,precomputed_user_path=None):
    global PRECOMPUTED_PICKLE
    global precomputed_ics

    print("Using Precalculated IC values")
    input_data = read_input(filepath)
    if precomputed_user_path is not None:
        PRECOMPUTED_PICKLE=precomputed_user_path
        precomputed_ics=unpickle_object(PRECOMPUTED_PICKLE)
    with open(GOA_UNIVERSE, "r") as fread:
        if not os.path.exists(UNIVERSE_COUNT_PICKLE):
            annotations = compile_universe(
                [x for x in fread.read().split("\n")])
            with open(UNIVERSE_COUNT_PICKLE, "wb") as fwrite:
                pickle.dump(annotations, fwrite)
        else:
            with open(UNIVERSE_COUNT_PICKLE, "rb") as fread:
                annotations = pickle.load(fread)
    print("Computing ICs")
    ic_data = compute_precalc(input_data, annotations, godag)
    print("Writing data file")
    dump_ic_data(
        ic_data, outputfile=filepath.split(
            "/")[-1].split(".")[0] + ".tsv"
    )
    # try:
    print("Plotting results")
    plot_density(ic_data, filepath.split("/")[-1].split(".")[0])

def precompute_data(annotpath,outputpath):
    global annotations
    if ".gaf" in annotpath:
        print("Parsing gaf file")
        if ".gz" in annotpath:
            with gzip.open(annotpath,"rt") as fread:
                annot_data=[x.split("\t")[4] for x in fread.read().split("\n") if x and (not x.startswith("!"))]        
        else:
            with open(annotpath,"r") as fread:
                annot_data=[x.split("\t")[4] for x in fread.read().split("\n") if x and (not x.startswith("!"))]        
    else:
        print("Reading '\n' separated GO list")
        with open(annotpath,"r") as fread:
            annot_data=[x for x in fread.read().split("\n") if x]

    input_data=process_file([[x,x] for x in set(annot_data)])
    
    annotations=compile_universe([x for x in annot_data if x])
    
    ic_data=compute_compute(input_data,annotations,godag)
    ic_dict={x[2]:[x[3],x[4],float(x[5])] for x in ic_data}
    print("Writing",outputpath)
    pickle_object(ic_dict,outputpath)

###############################################################################
# URLs
###############################################################################

# This may be deprecated some day
OBO_URL = "http://purl.obolibrary.org/obo/go.obo"
OBO_BASIC_URL = (
    "http://purl.obolibrary.org/obo/go-basic.obo"  # This may be deprecated some day
)

# Folder structure
DATA_FOLDER = new_folder("data/")
ANNOTATIONS_FOLDER = new_folder("annotations/")
PLOTS_FOLDER = new_folder("plots/")
TSVOUTPUTS_FOLDER = new_folder("outputs/")

# Data files
# We can use go-basic.obo or go.obo (go-basic is a simplified version)
if BASIC_OBO:
    OBO_PATH = get_remote_file(OBO_BASIC_URL, DATA_FOLDER + "go-basic.obo")
else:
    OBO_PATH = get_remote_file(OBO_URL, DATA_FOLDER + "go.obo")


DAG_PATH = DATA_FOLDER + "godag.pronto"
GOA_UNIVERSE=ANNOTATIONS_FOLDER+"goa_uniprot_all.universe"

# PICKLE
UNIVERSE_COUNT_PICKLE = DATA_FOLDER + "universe.pickle"
PRECOMPUTED_PICKLE = DATA_FOLDER+"goa_uniprot_all.pickle"

precomputed_ics = unpickle_object(PRECOMPUTED_PICKLE)

###############################################################################
# Parse OBO
###############################################################################
obo = parse_obo(OBO_PATH)

###############################################################################
# GO DAG
###############################################################################
if os.path.exists(DAG_PATH):
    godag = unpickle_object(DAG_PATH)
else:
    godag = Ontology(OBO_PATH)
    pickle_object(godag, DAG_PATH)
