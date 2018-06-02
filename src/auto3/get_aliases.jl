# parse the petsc_constants.jl file and get the typealiases
# this is portable (rather than printing their values, which is not)

# this requires petsc_constants to have a particular structure
# typealias name type
# global const thisname = val
# blank line here to separate from the next enum

"""
  Map (all possible) datatype strings to the printf format specifier
"""
global const FORMAT_SPEC = Dict{String, String}(
"String" => "%s",
"Cint" => "%i",
"Int32" => "%i",
"Uint32" => "%u",
"UInt32" => "%u",
)

"""
  Map enum name to datatype (string representation)
"""
function get_enum_types()
  f = open("../petsc_constants_gen.jl")

  enum_type_dict = Dict{String, String}()
  for line in eachline(f)
    if startswith(line, "typealias")
      words = split(line)
      @assert length(words) == 3
      enum_type_dict[ words[2] ] = words[3]
    end
  end

  close(f)
  return enum_type_dict
end

function get_enum_names()
  f = open("../petsc_constants_gen.jl")

  in_enum = false  # mode flag
  typealias_name = ""  # current enum name
  enum_name_dict = Dict{String, Vector{String}}()
  for line in eachline(f)
    # find the typealias, get all constant names following it (until space)
    if startswith(line, "typealias")
      in_enum = true
      words = split(line)
      typealias_name = words[2]
      
      if typealias_name in keys(enum_name_dict)
        error("repeat typealias named $typealias_name detected")
      end

      # make blank vector to populate later
      enum_name_dict[typealias_name] = String[]
    elseif in_enum && startswith(line, "global const")  # get an enum name
      words = split(line)
      push!(enum_name_dict[typealias_name], words[3])
    else  # end enum section
      in_enum = false
      typealias_name = ""
    end

  end  # end loop over lines

  close(f)
  return enum_name_dict
end


function check_enums(type_dict, name_dict)

  for i in keys(type_dict)
    @assert haskey(name_dict, i)
  end

  for i in keys(name_dict)
    @assert haskey(type_dict, i)
  end

  return nothing
end

"""
  This function writes the body of a C program that, when run, creates the
  petsc_constants_gen.jl file.
"""
function print_cbody(type_dict, name_dict, f::IO)

  check_enums(type_dict, name_dict)

  indent = "  "
  nl = "\\n"  # newline
  for i in keys(type_dict)
    datatype = type_dict[i]
    println(f, "\n", indent, "// typealias $i")
    println(f, indent, "fprintf(f, \"typealias $i $(datatype)$(nl)\");")

    # add quotes around strings
    if datatype == "String"
      mark = "\\\""
    else
      mark = ""
    end

    fmt = FORMAT_SPEC[datatype]
    for j in name_dict[i]
      println(f, indent, "fprintf(f, \"global const $j = $datatype( $mark$fmt$mark )$nl\", $j);")
    end  # end loop over enum names
    println(f, indent, "fprintf(f, \"\\n\");")
  end  # end loop over typealiases

  return nothing
end


function print_cheader(f::IO)

  indent = ""
  println(f, "#include <petscksp.h>\n")
  println(f, "int main(int argc, char **args)")
  println(f, "{")
  indent *= "  "
  println(f, indent, "PetscInitialize(&argc, &args, (char*)0, (char*)0);")
  println(f, "")
  println(f, indent, "FILE* f = fopen(\"petsc_constants_gen.jl\", \"w\");")
  println(f, "")

  return nothing
end

function print_cfooter(f::IO)
  
  indent = "  "
  println(f, "\n", indent, "PetscFinalize();")
  println(f, "\n",  indent, "return 0;")
  println(f, "}")

  return nothing
end


type_dict = get_enum_types()
name_dict = get_enum_names()

f = open("petsc_constants_gen.c", "w")
print_cheader(f)
print_cbody(type_dict, name_dict, f)
print_cfooter(f)
