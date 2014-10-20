import os
import os.path
import sys

# *****************************************************
#   check options and compatibility between modules
# *****************************************************

get_arch = 0
for x in sys.argv[1:]:

  if (x == "--get-arch"):
    print "Checking system architecture"
    get_arch = 1
    break

  if (x == "--with-chombo" or x == "--with-chombo:"): 
    print "Enabling Chombo support for AMR"
    for y in sys.argv[1:]:
      if (y == "--with-fd" or y == "--with-sb" or y == "--with-fargo"):
        print "! Incompatible modules, ",x," + ",y
        sys.exit(1)
    break

  elif (x == "--with-sb"): 
    print "Enabling support for shearing box module"
    for y in sys.argv[1:]:
      if (y == "--with-fd"):
        print "! Incompatible modules, ",x," + ",y
        sys.exit(1)
  
  elif (x == "--with-fd"): 
    print "Enabling support for finite difference module"

  elif (x == "--with-fargo"): 
    print "Enabling support for FARGO scheme"

  elif (x == "--no-curses"):
    print ""

  elif (x == "--help" or x == "-help"):
    print "Usage: python $sgPLUTO_DIR/setup.py [options]\n" 
    print "Here [options] can be:\n"
    print " --with-sb       enable the shearing box module."
    print " --with-fd       enable the finite difference module."
    print " --with-fargo    enable the FARGO-MHD module"
    print " --with-chombo   enable support for adaptive mesh refinement."
    print "                 (AMR) module using the Chombo library."
    print " --no-curses     disable ncurses library and use a"
    print "                 simpler text-based menu."
    sys.exit(1)

  else:
    print "! Unrecognized option '",x,"'"
    sys.exit(1)


print 'loading...'

# check if $sgPLUTO_DIR is defined

try:
  pluto_dir = os.environ["sgPLUTO_DIR"]
except:
  pluto_dir = os.getcwd()
#  print "---------------------------------------------------------------"
#  print " WARNING:"
#  print "  "
#  print " Your sgPLUTO_DIR shell variable does not seem to be defined. " 
#  print " Please set sgPLUTO_DIR as the directory name of your PLUTO    "
#  print " distribution before proceeding, e.g.\n"
#  print '    > setenv sgPLUTO_DIR "/home/user/PLUTO"\n'
#  print " if you're using tcsh, or \n"
#  print '    > export sgPLUTO_DIR="/home/user/PLUTO"\n'
#  print " if you're using bash."
#  print "---------------------------------------------------------------"
#  sys.exit()


print "PLUTO directory is defined as: ",pluto_dir
sys.path.append(pluto_dir+'/Tools/Python')
import ut
import configure

#sys.stdout.write('loading...')
#sys.stdout.flush()
#pluto_dir   = os.getcwd();
#pluto_dir   = os.path.abspath(pluto_dir)

configure.check(pluto_dir, get_arch)
if (get_arch): sys.exit(1)
#os.chdir(pluto_dir)
ut.main_menu(pluto_dir)

