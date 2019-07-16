#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jul 12 13:47:19 2019

@author: itodaiber
"""

file_input = "generated/limit.f90"
file_output = "compilable/limit_new.f90"

new_data = []

# add header
new_data = ["module limit_new\n",
            "  implicit none\n",
            "\n",
            "  private\n",
            "  \n",
            "  public :: get_associated_legendre_polynomial_limit\n",
            "\n",
            "  contains\n",
            "\n",
            "    subroutine io_error(str)\n",
            "      implicit none\n",
            "      ! dummy\n",
            "      character(len=*) :: str\n",
            "\n",
            "      print *, \"io ERROR\"\n",
            "      print *, \"  \", str\n",
            "    end subroutine\n",
            "\n",
            "    function cot(x)\n",
            "      implicit none\n",
            "      ! dummy\n",
            "      double precision :: x\n",
            "      double precision :: cot\n",
            "\n",
            "      cot = 1.0/tan(x)\n",
            "    end function\n",
            "\n",
            "    function csc(x)\n",
            "      implicit none\n",
            "      ! dummy\n",
            "      double precision :: x\n",
            "      double precision :: csc\n",
            "\n",
            "      csc = 1.0_8/sin(x)\n",
            "    end function\n",
            "\n",
            "    function get_associated_legendre_polynomial_limit(l,m) result(erg)\n",
            "      implicit none\n",
            "      ! dummy\n",
            "      integer(kind=4) :: l\n",
            "      integer(kind=4) :: m\n",
            "\n",
            "      double precision :: erg\n",
            "\n"]

def line_truncation(line, length):
    
    def add_space(line, n):
        line = (' ' * n) + line
        return line
    rv = []
    
    line_in = "".join(line)
    
    count = 0
    
    for i in line_in:
        if i == ' ':
            count += 1
        else:
            break
    
    buffer = []
    buffer_line = add_space("", count)
    
    print(len(line_in))
    
    if (len(line_in) > length):
        buffer = line_in.split()
        
        for b in buffer:            
            if (len(buffer_line) + len(b) < length):
                buffer_line = buffer_line + " " + "".join(b)
            else:
                if (len(buffer_line.lstrip()) > 0):
                    rv.append(buffer_line + " &\n")
                buffer_line = add_space("", count+4)
                buffer_line = buffer_line + "" + "".join(b)
            
        rv.append(buffer_line)
    else:
        rv = line
    return rv

with open(file_input, "r") as f:
    # read file line by line without '\n'-char
    data = f.read().splitlines()
     
    flag_at_erg_line = False
    buffer = []
    
    for line in data:
        # adapt some function calls
        line = line.replace("\\", "")   # remove backslash
        line = line.replace("endselect","end select")
        line = line.replace("*********************", "!*********************")
        
        # add indentation
        line = "      " + line
        
        # find "erg="-lines and append all further lines up to "case"-line
        if (line.find("erg=") > 0):
            flag_at_erg_line = True
            buffer.append(line)
        elif(line.find("case") > 0):
            if (flag_at_erg_line is True):
                flag_at_erg_line = False
                b = line_truncation(buffer, 80)
                new_data.append(" ".join(b) + "\n")
                buffer=[]
                new_data.append(line + "\n")
            else:
                new_data.append(line + "\n")
        else:
            if (flag_at_erg_line is True):
                buffer.append(line.lstrip())
            else:
                new_data.append(line + "\n")

# add footer
new_data.append("".join(["  end function get_associated_legendre_polynomial_limit\n",
                         "end module\n"]))

# write to file
with open(file_output, "w") as f:
    f.writelines(new_data)

print("Ready")