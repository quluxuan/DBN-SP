import math
def pprint(file,obj,level=0,indent="  "):
    """Pretty print a dictionary
    """
    try:
        items=obj.items()
    except AttributeError:
        file.write(repr(obj))
    else:
        writeln_indent(file,"{",level*indent)
        for k,val in sorted(items):
            write_indent(file,repr(k)+" : ",(level+1)*indent)
            pprint(file,val,level+1)
            writeln_indent(file,",","")
        write_indent(file,"}",level*indent)
def safe_exponent(sc):
    """Gives the safe (against overflows) exponent of sc (base 2).

    In case of overflow, give the approximate value.
    """
    try:
        prob=str(2**sc)
    except OverflowError:
        sc10=math.log(2,10)*sc
        prob=str(10**(sc10-int(sc10)))+"e+"+str(int(sc10))
    return prob
def writeln_indent(file,txt,ind):
    """write an indented string followed by EOL
    """
    write_indent(file,txt,ind,ln=True)
def write_indent(file,txt,ind,ln=False):
    """write a string to a file with indent
    """
    file.write(ind)
    file.write(txt)
    if ln:
        file.write("\n")