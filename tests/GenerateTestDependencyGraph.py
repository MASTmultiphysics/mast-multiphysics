#!/usr/bin/env python

import os
import sys

cmakeFilePath = "/home/neiferd/OpenSourceProjects/MAST-Library/tests/CMakeLists.txt"


class CMakeStructure(object):
    def __init__(self, root_cmake_lists_path):
        self.path = root_cmake_lists_path
        self.lines = self.read_cmake(root_cmake_lists_path)
        self.subdirectories = []
        self.tests = {}

    def get_dependency_structure(self):
        self.get_subdirectories()

        for subdir in self.subdirectories:
            lines = self.read_cmake(subdir)
            a.get_tests(lines)
            a.get_tests_properties(lines)

    def write_graph_code(self):
        # outpath = os.path.dirname(cmakeFilePath)+"/auto_test_dependency_graph.dot"
        outpath = "./test_dependency_graph.dot"
        with open(outpath, 'w') as graph:
            graph.write("digraph test_dependencies\n{\n")
            for test_name,test in self.tests.items():
                for fr in test.fixtures_required:
                    graph.write('    "%s" -> "%s"\n'%(fr,test_name))
            graph.write("}")
        return outpath

    def read_cmake(self, file_path):
        with open(file_path) as file:
            lines = file.readlines()
        lines = [x.strip() for x in lines if ((not x.startswith("#")) and (len(x.strip()) > 0))]
        return lines

    def get_subdirectories(self, file_path=None):
        if file_path is None:
            file_path = self.path
        lines = self.read_cmake(file_path)
        for line in lines:
            if line.startswith("add_subdirectory("):
                subdir_path = os.path.dirname(file_path) + "/" + line[line.find("(") + 1:line.find(")")] + "/CMakeLists.txt"
                self.subdirectories.append(subdir_path)
                self.get_subdirectories(subdir_path)

    def get_tests(self, lines):
        same_test = False
        for line in lines:
            if not same_test:
                if line.startswith("add_test("):
                    test = Test()

                    data = line.lstrip("add_test(").split()
                    for i in range(0, len(data), 2):
                        if data[i].lower() == "name":
                            test.name = data[i+1].rstrip(")")
                        elif data[i].lower() == "command":
                            test.command = data[i+1].rstrip(")")

                    if line.endswith(")"):
                        same_test = False
                        self.tests[test.name] = test
                    else:
                        same_test = True

            else:
                data = line.lstrip("add_test(").split()
                for i in range(0, len(data), 2):
                    if data[i].lower() == "name":
                        test.name = data[i + 1]
                    elif data[i].lower() == "command":
                        test.command = data[i + 1]
                if line.endswith(")"):
                    same_test = False
                    self.tests[test.name] = test
                else:
                    same_test = True

    def get_tests_properties(self, lines):
        same_test = False
        for line in lines:
            if not same_test:
                if line.startswith("set_tests_properties("):
                    data = line.lstrip("set_tests_properties(").split()
                    test_name = data[0]
                    for i in range(1,len(data)):
                        if data[i].lower() == "fixtures_required":
                            if data[i+1].startswith('"') and data[i+1].endswith('"'):
                                frdata = data[i+1][1:-1].split(";")
                                for val in frdata:
                                    self.tests[test_name].fixtures_required.append(val)
                            else:
                                self.tests[test_name].fixtures_required.append(data[i+1].rstrip(")"))
                        elif data[i].lower() == "fixtures_setup":
                            self.tests[test_name].fixtures_setup.append(data[i + 1].rstrip(")"))
                        elif data[i].lower() == "fixtures_cleanup":
                            self.tests[test_name].fixtures_cleanup.append(data[i + 1].rstrip(")"))
                    if line.endswith(")"):
                        same_test = False
                    else:
                        same_test = True
            else:
                data = line.lstrip("set_tests_properties(").split()
                for i in range(0, len(data)):
                    if data[i].lower() == "fixtures_required":
                        if data[i + 1].startswith('"') and data[i + 1].endswith('"'):
                            frdata = data[i + 1][1:-1].split(";")
                            for val in frdata:
                                self.tests[test_name].fixtures_required.append(val)
                        else:
                            self.tests[test_name].fixtures_required.append(data[i + 1].rstrip(")"))
                    elif data[i].lower() == "fixtures_setup":
                        self.tests[test_name].fixtures_setup.append(data[i + 1].rstrip(")"))
                    elif data[i].lower() == "fixtures_cleanup":
                        self.tests[test_name].fixtures_cleanup.append(data[i + 1].rstrip(")"))
                if line.endswith(")"):
                    same_test = False
                else:
                    same_test = True


class Test(object):
    def __init__(self):
        self.name = None
        self.command = None
        self.fixtures_setup = []
        self.fixtures_cleanup = []
        self.fixtures_required = []
        self.labels = []
        self.processors = None


if __name__ == "__main__":
    if len(sys.argv)>=2:
        cmakeFilePath = sys.argv[1]
    else:
        cmakeFilePath = "./CMakeLists.txt"
    a = CMakeStructure(cmakeFilePath)
    a.get_dependency_structure()
    outpath = a.write_graph_code()

    try:
        os.system("dot -Tpng %s -o %s"%(outpath, outpath[:-3]+"png"))
        os.system("dot -Teps %s -o %s" % (outpath, outpath[:-3] + "eps"))
    except:
        print("Warning: Error during conversion from dot to png and/or eps.  Is graphviz installed?")
