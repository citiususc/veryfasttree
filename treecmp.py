#! /usr/bin/env python
import sys
try:
	from ete3 import Tree
except ImportError:
	print("Module ete3 required. Please install with pip install ete3 or easyinstall ete3", file=sys.stderr)
	exit(-1)

if len(sys.argv) <= 2:
	exit("use " + sys.argv[0] + " tree1 tree2")

try:
	tree1 = Tree(sys.argv[1])
	tree2 = Tree(sys.argv[2])

	rf, max_parts = tree1.robinson_foulds(tree2, unrooted_trees=True)[:2]
	accuracy = round((1 - (rf / max_parts)) * 100, 2)
	print("Accuracy: %.2f%%" % accuracy)

except Exception as e:
	print(str(e), file=sys.stderr)
