#!/usr/bin/env python

import os
import sys

colormap = {"Passed"    :   "green",
            "Not Run"   :   "yellow",
            "Failed"    :   "red"}

status = ["Passed", "Not Run", "Failed"]

def generate_result_graph(dotPath, resultsPath):
    # Get Test Results
    with open(resultsPath, 'r') as resf:
        resLines = resf.readlines()

    tests = {}
    for i in range(len(resLines)):
        line = resLines[i].strip()
        if ((" Test " in line) and ("#" in line)):
            data = line.split()
            test_name = data[3]
            for s in status:
                if s in line:
                    result = s
                    break
            tests[test_name] = result


    # Generate color code lines for dot graph
    with open(dotPath, 'r') as dotf:
        dotLines = dotf.readlines()

    cdotLines = [dotLines[0], dotLines[1], "    node [style = filled];\n"]
    for test_name, result in tests.items():
        cdotLines.append('    "%s" [color=%s];\n'%(test_name, colormap[result]))
    cdotLines += dotLines[2:]


    # Write new dot to file
    # outPath = os.path.dirname(dotPath)+"/results_"+os.path.basename(dotPath)
    outPath = "./test_results_graph.dot"
    with open(outPath, 'w') as f:
        [f.write(line) for line in cdotLines]

    try:
        os.system("dot -Tpng %s -o %s"%(outPath, outPath[:-3]+"png"))
        os.system("dot -Teps %s -o %s" % (outPath, outPath[:-3] + "eps"))
    except:
        print("Warning: Error during conversion from dot to png and/or eps.  Is graphviz installed?")


if __name__ == "__main__":
    print("\nGENERATING RESULTS GRAPH\n")
    dotPath = sys.argv[1]
    resultsPath = sys.argv[2]

    generate_result_graph(dotPath, resultsPath)