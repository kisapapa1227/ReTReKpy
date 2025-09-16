###### セル１
import utils
utils.get_default_config()
import sys
import json

from rdkit import Chem
from rdkit.Chem import Draw
from PIL import Image, ImageDraw
import os

####### セル２（分子の入力）
from utils import *
import os
from reimplemented_libraries import CxnUtils
from model_modules import load_model, predict_templates

config = get_default_config()
# Create save directory
name_stem = config["target"].split('/')[-1].split('.')[0]
os.makedirs(config["save_result_dir"], exist_ok=True)

_image_dir="/var/www/html/public/images"

# Setup logger
#level = DEBUG if args.debug else INFO
#logger = get_logger(level, config["save_result_dir"])

config["rollout_model"]="model/variables"
config["rollout_rules"]="data/reaction_template_uspto_filtered.sma"
config["expansion_model"]="model/variables"
config["expansion_rules"]="data/reaction_template_uspto_filtered.sma"
gateway = CxnUtils(config['rollout_rules'])

# data preparation
smiles = sys.argv[1]
#smiles = 'CCCOCO'
target_mol = Chem.MolFromSmiles(smiles)

expansion_rules = ReactionUtils.get_reactions(config['expansion_rules'], config['save_result_dir'])
rollout_rules = ReactionUtils.get_reactions(config['rollout_rules'], config['save_result_dir'])
with open(config['starting_material'], 'r') as f:
    start_materials = set([s.strip() for s in f.readlines()])

#with open(config['intermediate_material'], 'r') as f:
##    intermediate_materials = set([s.strip() for s in f.readlines()])

intermediate_materials = set()

expansion_model = load_model('expansion', config, class_num=len(expansion_rules))
rollout_model = load_model('rollout', config, class_num=len(rollout_rules))

#template_scores = json.load(open(config["template_scores"], "r"))
template_scores = {}
in_scope_model = load_model('in_scope', config)

####### セル３
start_materials_=list(start_materials)
start_materials_
####### セル5(ここが詳細設定の入力欄)
# 引数がない場合はデフォルトのリストを使用
default_weights = [1.0, 1.0, 1.0, 1.0, 1.0, 1.0]
if len(sys.argv) > 3:
    json_weights = sys.argv[3]
    try:
        # JSON形式の文字列をPythonのリストに変換
        knowledge_weights = json.loads(json_weights)
    except json.JSONDecodeError:
        # JSON形式として解釈できない場合はデフォルト値
        knowledge_weights = default_weights
else:
    knowledge_weights = default_weights

save_tree = sys.argv[4] == 'True' if len(sys.argv) > 4 else False
expansion_num = int(sys.argv[5]) if len(sys.argv) > 5 else 50
cum_prob_mod = sys.argv[6] == 'True' if len(sys.argv) > 6 else False
chem_axon = sys.argv[7] == 'True' if len(sys.argv) > 7 else False
selection_constant = int(sys.argv[8]) if len(sys.argv) > 8 else 10
time_limit = int(sys.argv[9]) if len(sys.argv) > 9 else 0
csrf_token = sys.argv[10]
csv = sys.argv[11]
substance = sys.argv[12]
uid = sys.argv[13]
upper_time = int(sys.argv[14])

config['knowledge']=['all']
config['knowledge_weights']=knowledge_weights
config['save_tree']=save_tree
config['expansion_num']=expansion_num
config['cum_prob_mod']=cum_prob_mod
config['chem_axon']=chem_axon
config['selection_constant']=selection_constant
config['time_limit']=time_limit
config['debug']=True

####### セル６(1)
from mcts_main import Mcts
from mcts_modules import print_route_HTML, save_route
from myMcts_modules import myPrint_route_HTML

import logging
logging.basicConfig(
    stream=sys.stdout,
    level=logging.DEBUG,
    force=True)

logger = logging.getLogger()
# main process

mcts = Mcts(target_mol, expansion_rules, rollout_rules, start_materials, intermediate_materials, template_scores, config)

#logger.info(f"[INFO] knowledge type: {config['knowledge']}")
#logger.info("[INFO] start search")

import time
import subprocess

# logging設定: DEBUG メッセージを非表示に
logging.basicConfig(level=logging.INFO, format='%(levelname)s:%(message)s')
# PILのlogging levelをWARNINGに設定する
logging.getLogger('PIL').setLevel(logging.WARNING)
logger = logging.getLogger(__name__)
route_num = int(sys.argv[2]) if len(sys.argv) > 2 else 3
#komai

out_dir=f"{_image_dir}/{uid}"
if len(substance)<1:# when the name is omitted
    dt_now=datetime.datetime.now()
    substance="search_"+str(dt_now.year)+str(dt_now.month).zfill(2)+str(dt_now.day).zfill(2)+str(dt_now.hour).zfill(2)+str(dt_now.minute).zfill(2)

fn=out_dir+"/"+substance+".txt"

if not os.path.isdir(out_dir):
     os.mkdir(out_dir)

if os.path.isfile(fn):
    os.remove(fn)

#    sleep(1)
while os.path.isfile(fn):
    sleep(1)

with open(fn,"w") as fp:
    fp.write(f"{smiles}\n")

init = time.time()

route=0
while True:
#for route in range(route_num):
    now = time.time()
    if upper_time!=0 and (now-init)>upper_time*60:
        break
    if route_num!=0 and route == route_num:
        break
    route = route + 1
    logger.setLevel(logging.WARNING)
    start = time.time()
    leaf_node, is_proven = mcts.search(expansion_model, rollout_model, in_scope_model, logger, gateway=gateway, time_limit=config['time_limit'])
    elapsed_time = time.time() - start
    logger.info(f"[INFO] done in {elapsed_time:5f} s")
    logger.setLevel(logging.INFO) 

    nodes = []
    while leaf_node.parent_node is not None:
        nodes.append(leaf_node)
        leaf_node = leaf_node.parent_node
    else:
        nodes.append(leaf_node)

    myPrint_route_HTML(nodes, is_proven, logger, route_num, smiles, route, json_weights, save_tree, expansion_num, cum_prob_mod, chem_axon, substance, selection_constant, time_limit, csrf_token, fn, image_dir=out_dir);

db=f'tmp{uid}.db'

if os.path.isfile(db):
    os.remove(db)

#tt=['python3','make_reports/addDb.py','-u',uid,'-d',_image_dir,'-database',db,"-n",str(route)]
#std=subprocess.run(tt,capture_output=True,text=True).stdout

#with open(fn+"_log","a") as fp:
#    fp.write("-----> addDb"+str(route)+"\n");
#    for t in std.split("\n"):
#        fp.write("ok:"+t+"\n")

#tt=['python3','make_reports/readDb.py','-u',uid,'-id','1','-d',_image_dir+'/'+uid+'/report/','-database',db,"-chem",substance,"-force"]

#std=subprocess.run(tt,capture_output=True,text=True).stdout
#with open(fn+"_log","a") as fp:
#    fp.write("-----> readDbDb\n");
#    for t in std.split("\n"):
#        fp.write("ok:"+t+"\n")

with open(fn,"a") as fp:
    fp.write("searched:"+substance+"\n")

#time.sleep(1);
#with open(fn+"_log","w") as fp:
#with open(fn+"_log","a") as fp:
#    fp.write("routeAnalyzer\n")
#    fp.write(std+"\n")
