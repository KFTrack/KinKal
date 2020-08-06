#!/usr/bin/env python
#
# Define some helper functions and one helper class, all used to build the shared libraries.
#
# The functions are:
# def validateEnvironment():
#     - Check that some important environment variables are set.  Exit if they are not.
#
# def validateDebugLevel( ):
#     - Check that the debug level has one of the supported values
#
# def defineMergeFlags( debug_level):
#     - Define the compiler and linker flags; some flags depend on the debug leve.
#
# def defineExportedOSEnvironment():
#     - Prepare a set a of environment variables, plus their values, that will
#       be added to the scons environment.
#
# def locateSConscriptFiles(sourceRoot):
#     - Walk the directory tree rooted at sourceRoot and return a list of the paths to
#       all SConscript files under that root.
#
# def dispatchSConscriptFiles( ss, sourceRoot, build_base ):
#     - Given a list of SConscript files, ss, tell scons to build the targets specified
#       in this files and to put the outputs into the defined spot in the build heirarchy.
#
# The class build_helper has four user callable entry points:
#
#  1) make_mainlib( self, userlibs, cppf=[], pf=[], addfortran=False ):
#     Find all non-plugin .cc files in the current directory, compile
#     them and add them to a library named after the path to the current diretory.
#     The user must provide a link list (userlibs), which may be empty.
#     Optionally provide:
#        - additional compiler flags (cppf)
#        - additional flags to the parse flags feature
#     Optionally compile all fortran files and add them to the library.
#
#  2) make_plugin( self, cc, userlibs, cppf = [], pf = []):
#     Make single plugin library, specified by the named .cc file.
#     The user must provide a link list (userlibs), which may be empty.
#     The library is named after that path to the source file.
#     The optional arguments are the same as for 1)
#
#  3) make_plugins( self, userlibs, exclude_cc = [], cppf = [], pf = [] ):
#     Find all plugin files in the current directory and build them against the
#     same link list; the link list may be empty.  Exclude any plugins specified
#     on the exclude_cc list.  For each plugin this makes a shared library
#     named after the path to the source file.
#     This does not make the _dict and _map plugins.
#


from glob import glob
import os, re, string
from timeit import default_timer as timer
import sys
import subprocess
from pathlib import Path
# Check that some of the required environment variables have been set.
def validateEnvironment():
    if not 'PACKAGE_SOURCE' in os.environ:
        sys.exit('You have not specified PACKAGE_SOURCE for this build.\nExiting.')
    if not 'DEBUG_LEVEL' in os.environ:
        sys.exit('You have not specified DEBUG_LEVEL for this build.\nExiting.')
    if not 'BUILD_BASE' in os.environ:
        sys.exit('You have not specified BUILD_BASE for this build.\nExiting.')

# Check that the debug level is one of the known values
def validateDebugLevel():
    level = os.environ['DEBUG_LEVEL']
    allowedLevels = ['prof', 'debug' ]
    if not level in allowedLevels:
        print ('Unrecognized value for --debuglevel ', level)
        print ('   The value must be one of the allowed levels: ', str(allowedLevels) )
        raise Exception('Illegal value for --debuglevel')
        pass
    print ("Recognized debug level is: ", level)
    return level

# Define the compiler and linker options.  Scons knows how to sort out which is which.
# These are given to scons using its Evironment.MergeFlags call.
def defineMergeFlags( debug_level):
#   commonFlags = [ '-std=c++17', '-rdynamic', '-Wall', '-Wno-unused-local-typedefs', '-g', '-Werror',
   commonFlags = [ '-std=c++17', '-Wall', '-Wno-unused-local-typedefs', '-g', '-Werror',
                   '-gdwarf-2', '-Werror=return-type', '-Winit-self', '-Woverloaded-virtual', '-fPIC', '-ftrapping-math' ]
                   # Linker flags
                   #'-Wl,--no-undefined',
#                   '-Wl,--as-needed' ]

   # Hopefully the need for this will go away in later compiler versions.
   gcc_version = os.environ['GCC_VERSION']
   if gcc_version == 'v7_3_0':
       commonFlags = [ commonFlags, '-Wno-strict-overflow' ]
       pass

   extraFlags = []
   if debug_level == 'prof':
       extraFlags = [ '-O3', '-fno-omit-frame-pointer', '-DNDEBUG' ]
       pass
   if debug_level == 'debug':
       extraFlags = [ '-O0'  ]
#       extraFlags = [ '-O0', '-ferror-limit=0'  ] sadly this isn't supported by gcc
       pass
   return [ commonFlags, extraFlags ]

# Prepare some shell environment in a form to be pushed into the scons environment.
def defineExportedOSEnvironment():
    osenv = {}
    for var in [ 'LD_LIBRARY_PATH',  'GCC_FQ_DIR',  'PATH', 'PYTHONPATH',  'ROOTSYS', 'PYTHON_ROOT', 'PYTHON_DIR' ]:
        if var in os.environ.keys():
            osenv[var] = os.environ[var]
            pass
        pass
    return osenv

# Walk the directory tree to locate all SConscript files.
def locateSConscriptFiles(sourceRoot):
    ss=[]
    exclude = [Path(os.environ['BUILD_BASE'])]

    for root, _, files in os.walk(sourceRoot):
        candidate_path = os.path.join(root,'SConscript')

        if 'SConstruct' in files and root != sourceRoot:
            # exclude any SConscripts found under a build directory
            rootp = Path(root)
            if rootp not in exclude:
                exclude.append(rootp)
            continue

        if 'SConscript' in files:
            ss.append(candidate_path)

    # remove SConscript paths found under any build directories
    sconscripts = []
    for cand_path in ss:
        append = True
        for ex in exclude:
            if ex in Path(cand_path).parents:
                append = False
        if append:
            sconscripts.append(cand_path)
    return sconscripts

# Tell scons to do the work found in the set of SConscript files.
def dispatchSConscriptFiles( env, ss, sourceRoot, build_base ):
    for sourceFile in ss:
        relpath = os.path.relpath(sourceFile, sourceRoot )
        tokens = relpath.split('/')
        tokens.pop()
        sep = '/'
        objPath = build_base + '/' + sep.join(tokens)
        env.SConscript ( sourceFile, variant_dir=objPath,)

# An instance of this class available to the SConscript files via the scons environment.
class build_helper:
    """build_helper: class to build libraries from source"""
#   This behaves like c++ static data member and is initialized at construction time.
    def __init__(self,env):
        self.env       = env
        self.buildBase = os.environ['BUILD_BASE']
#
#   Accesor
#
    def base(self):
        return self.buildBase
#
#   Build the name of the shared library into which non-plugin compiled code will be inserted.
#   Two versions: with and without the '#/lib' path prefix.
#
    def libname(self):
        relpath = os.path.relpath('.',self.buildBase)
        tokens = relpath.split('/')
        if len(tokens) > 1:
            if tokens[len(tokens)-1] == 'src':
                tokens.pop()
                pass
            pass
        sep = '_'
        return sep.join(tokens)

    def prefixed_libname(self):
        return '#/lib/' + self.libname()
#
#   Build the name of the shared library into which plugin code will be inserted.
#   Two versions: with and without the '#/lib' path prefix.
#
    def plugin_libname(self,sourcename):
        return self.libname() + '_' + sourcename[:sourcename.find('.cc')]
    def prefixed_plugin_libname(self,sourcename):
        return '#/lib/' + self.plugin_libname(sourcename)
#
#   Build a list of plugins to be built.
#
    def plugin_cc(self):
        return ( self.env.Glob('*_module.cc',  strings=True) +
                 self.env.Glob('*_service.cc', strings=True) +
                 self.env.Glob('*_source.cc',  strings=True) )

#
#   Build a list of .cc files that are not plugings; these go into the
#   library named after the directory.
#
    def non_plugin_cc(self):
        tmp = non_plugin_cc = self.env.Glob('*.cc', strings=True)
        for cc in self.plugin_cc(): tmp.remove(cc)
        for cc in self.unittest_cc(): tmp.remove(cc)
        return tmp
#
#   Names need to build the _dict and _map libraries.
#
    def dict_tmp_name(self):
        relpath = os.path.relpath('.',self.buildBase)
        return '#/tmp/' + relpath + '/' + self.libname() + '_dict.cpp'

    def map_tmp_name(self):
        relpath = os.path.relpath('.',self.buildBase)
        return '#/tmp/' + relpath + '/' + self.libname() + '_map.cpp'

    def dict_libname(self):
        relpath = os.path.relpath('.',self.buildBase)
        return self.libname() + '_dict'

    def map_libname(self):
        relpath = os.path.relpath('.',self.buildBase)
        return self.libname() + '_map'

    def prefixed_dict_libname(self):
        return '#/lib/' + self.dict_libname()

    def prefixed_map_libname(self):
        return '#/lib/' + self.map_libname()
#
#   Make the main library.
#
    def make_mainlib( self, userlibs, cppf=[], pf=[], addfortran=False ):
        non_plugin_cc = self.non_plugin_cc()
        if addfortran:
            fortran = self.env.Glob('*.f')
            non_plugin_cc = [ non_plugin_cc, fortran]
            pass
        libs = []
        if non_plugin_cc:
            libs = self.env.SharedLibrary( self.prefixed_libname(),
                                    non_plugin_cc,
                                    LIBS=[ userlibs ],
                                    CPPFLAGS=cppf,
                                    parse_flags=pf
                                )
            #libs = [ self.libname() ]
            pass
        return libs
#
#   Make one plugin library ( but does not work for _dict and _map plugins )
#
    def make_plugin( self, cc, userlibs, cppf = [], pf = []):
        self.env.SharedLibrary( self.prefixed_plugin_libname(cc),
                                cc,
                                LIBS=[ userlibs, ],
                                CPPFLAGS=cppf,
                                parse_flags=pf
                            )
#
#   Make all plugin libraries, excluding _dict and _map; this works if
#   all libraries need the same link list.
#
    def make_plugins( self, userlibs, exclude_cc = [], cppf = [], pf = [] ):
        plugin_cc = self.plugin_cc()
        for cc in exclude_cc: plugin_cc.remove(cc)
        for cc in plugin_cc:
            self.env.SharedLibrary( self.prefixed_plugin_libname(cc),
                                    cc,
                                    LIBS=[ userlibs ],
                                    CPPFLAGS=cppf,
                                    parse_flags=pf
            )

#
#   Make the dictionary and map plugins.
#
    def make_dict( self ):
        cmd=self.env.Command('Dict.cc',['KKHitInfo.hh','KKBFieldInfo.hh','KKMatInfo.hh','KTrajInfo.hh','LinkDef.h'], 
            'rootcling -f {build}/UnitTests/Dict.cc {src}/UnitTests/KKHitInfo.hh {src}/UnitTests/KKMatInfo.hh {src}/UnitTests/KTrajInfo.hh {src}/UnitTests/KKBFieldInfo.hh {src}/UnitTests/LinkDef.h && mv {build}/UnitTests/Dict_rdict.pcm {build}/lib/Dict_rdict.pcm'.format(
                src=os.environ['PACKAGE_SOURCE'], build=os.environ['BUILD_BASE']
            ))
        #self.env.Depends('UnitTests', 'Dict.cc')
        return cmd
#
#   Build a list of the source filenames for unit tests to be built and run
#
    def unittest_cc(self):
        return ( self.env.Glob('*_unit.cc',  strings=True) )

#
#   filename of executable to be made from the unit test source.
#   Just the filename without other path elements.
#
    def executable_name(self, name ):
       return name[:name.find('_unit.cc')]

#   Build all unit test
#
    def make_unit_tests( self, linkLists ):
        unitTests = self.unittest_cc()

        test_alias = self.env.Alias('test', '', self.run_unit_tests())
        self.env.AlwaysBuild(test_alias)
        
        pgs = []
        for u in unitTests:
            exe = self.executable_name(u)
            libs  = linkLists.get(exe,[])
            pg = self.env.Program(
                target = "#/bin/"+exe,
                source = [ u ],
                LIBS   = libs
            )
            self.env.Depends(test_alias, pg)
            pgs.append(pg[0])
        return pgs
#
#   Run all unit tests.
#

    def run_unit_tests(self):
        unitTests = self.unittest_cc()

        def run_unit_test( executable, print_output=True):
            cmd = "./bin/" + executable

            p=subprocess.Popen(cmd, stdout=subprocess.PIPE, env={**os.environ}, shell=True)
            (output, err) = p.communicate()
            p_status = p.wait()

            if print_output and output is not None:
                print (output.decode())
            if print_output and err is not None:
                print (err.decode())

            return p_status, output, err


        def runner(*args,**kwargs):
            statuses = []
            tests_failing = 0
            for u in unitTests:
                exe = self.executable_name(u)
                print( )
                print("Running %s..." % exe)

                start = timer()
                rc, _, _ = run_unit_test( exe )
                end = timer()

                status = ' %.2f s [   OK   ]' % (end-start)
                if rc != 0:
                    tests_failing += 1
                    status = ' %.2f s [ FAILED ]' % (end-start)
                
                statuses.append('- '+exe + ' ' * (30 - len('- '+exe)) + status)

                #print ( "   Done: ", exe, "  status: ", status)
            print( )
            print( )
            print ("UNIT TEST SUMMARY")
            print ('-' * (30 + len(' %.2f s [ ...... ]')))
            for l in statuses:
                print (l)
            print( )
            print ( '%2d / %2d tests passed.' % (len(unitTests)-tests_failing, len(unitTests)))
            print( )
            print( )

            return tests_failing

        return runner
