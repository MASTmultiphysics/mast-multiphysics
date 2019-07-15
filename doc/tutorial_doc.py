import os, sys


def sanitize_comment(l):
    l = l.strip()

    if (len(l) > 1 and l[0:2] ==  "//"):
        l = " * " + l.strip("//") + " \n"

    return l

def if_line_is_comment(l):
    ln = l.strip()
    if ((len(ln) > 1 and ln[0:1] == "* ") or (len(ln) > 1 and ln[0:2] == "//")):
        return True
    else:
        return False


def parse_example(physics, d, physics_dir, doc_dir):
    # the examples are structured such that each directory has one main file named example_XX.cpp,
    # where, XX is the number of the example.
    source_abs_path = os.path.join( physics_dir, d, d+".cpp")
    output_abs_path = os.path.join(     doc_dir, physics+"_"+d+".dox")
    assert(os.path.isfile(source_abs_path))

    # now read the file and translate the documentation
    example_source_file = open(source_abs_path, "r")
    example_output_file = open(output_abs_path, "w")
    translate           = False
    code_mode           = False


    for l in example_source_file:
        if "BEGIN_TRANSLATE" in l:
            ln = l.split("BEGIN_TRANSLATE")
            if len(ln) > 1:
                example_output_file.write("/*! \n * \page " + physics+"_"+d + ln[1] + " \n")
            else:
                example_output_file.write("/*! \n * \page " + physics+"_"+d + " \n")
            translate = True
        elif "END_TRANSLATE" in l:
            if code_mode:
                example_output_file.write(" * \endcode \n")
            example_output_file.write("*/ \n")
            translate = False
        elif translate:
            if if_line_is_comment(l):
                if code_mode:
                    example_output_file.write(" * \endcode \n")
                    code_mode = False
                example_output_file.write(sanitize_comment(l))
            else:
                if not code_mode:
                    example_output_file.write(" * \code{.cpp} \n")
                    code_mode = True
                if (len(l) > 0):
                    example_output_file.write(l)

    example_source_file.close()
    example_output_file.close()



if __name__ == "__main__":

    # the call should be something like $PYTHON tutorial_doc.py $SOURCE_DIR
    assert(len(sys.argv) == 3)

    # terminal output to show during CMake execution.
    print(sys.argv[0])
    print("Processing examples/*/example_XX/* files.")

    examples_dir = sys.argv[1]
    doc_dir      = sys.argv[2]

    assert(os.path.isdir(examples_dir))
    assert(os.path.isdir(doc_dir))

    # make sure the subdirectories exists
    assert(os.path.isdir(examples_dir))
    assert(os.path.isdir(doc_dir))

    # get the list of example directories to be processed
    for physics in os.listdir(examples_dir):
        if (os.path.isdir(os.path.join(examples_dir, physics))):
            physics_dir = os.path.join(examples_dir, physics)
            for ex in os.listdir(physics_dir):
                ex_dir = os.path.join(physics_dir, ex)
                if ("example_" in ex and os.path.isdir(ex_dir)):
                    parse_example(physics, ex, physics_dir, doc_dir)

