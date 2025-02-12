
import os
import math
from SCons.Script import *
from SCons.Node.FS import File
from SCons import Action

"""
Scons Unity Build Generator

Provides several generators for SCons to combine multiple source files into a bigger
one to reduce compilation time, so called "unity builds". This is achieved by generating
unity source files which in term include the actual source files and compile them using
one of the existing SCons builders.

Usage
-----
In order to use this, just place it inside your `site_scons/site_tools` folder, enable it by
adding "unity_build" to the tools when constructing your Environment and replace invocations
of the Program/Library/SharedLibrary/StaticLibrary builders with their Unity... counterpart:

env = Environment(tools = ['default', 'unity_build'])

source_files = ...

env.UnityProgram(
    target = 'my_program',
    source = source_files,
    ...
)

The tool will generate an amount of unity source files and invoke the Program builder on these,
forwarding any other arguments you passed.

Other Options
------------
You can control the behaviour of the builder using several Environment options:
env['UNITY_CACHE_DIR'] = '.unity' # Directory where the unity sources are stored.
                                  # can be either a string or a Dir() node.
env['UNITY_MAX_SOURCES'] = 15     # Maximum number of source files per unity file.
env['UNITY_MIN_FILES'] = env.GetOption('num_jobs')
                                  # Minimum number of unity files to generate (if possible).
                                  # Defaults to the number of jobs passed to SCons.
env['UNITY_DISABLE'] = False      # Set to True to completely disable unity builds. The commands
                                  # will simply pass through their options to the regular builders.

Additionally any generator can be passed a `cache_dir` to overwrite the value from the Environment.
"""

def exists(env : Environment):
    return True

def generate(env : Environment):
    env.AddMethod(_make_generator(env.Program), 'UnityProgram')
    env.AddMethod(_make_generator(env.Library), 'UnityLibrary')
    env.AddMethod(_make_generator(env.StaticLibrary), 'UnityStaticLibrary')
    env.AddMethod(_make_generator(env.SharedLibrary), 'UnitySharedLibrary')

    # build for generating the unity source files
    unity_source_builder = env.Builder(
        action = Action.Action(_generate_unity_file, _generate_unity_file_msg)
    )
    env.Append(BUILDERS = {'UnitySource': unity_source_builder})

    env.SetDefault(UNITY_CACHE_DIR = '.unity')
    env.SetDefault(UNITY_MAX_SOURCES = 15)
    env.SetDefault(UNITY_MIN_FILES = env.GetOption('num_jobs'))
    env.SetDefault(UNITY_DISABLE = False)

def _make_generator(base_generator):
    def generator(env, source, target, cache_dir = None, *args, **kwargs):
        if env['UNITY_DISABLE']:
            return base_generator(target = target, source = source, *args, **kwargs)
        unity_source_files = []
        source_files, other_nodes = _flatten_source(source)

        max_sources_per_file = max(1, math.ceil(len(source_files) / env['UNITY_MIN_FILES']))
        sources_per_file = min(max_sources_per_file, env['UNITY_MAX_SOURCES'])
        
        num_unity_files = math.ceil(len(source_files) / sources_per_file)

        if not cache_dir:
            cache_dir = env['UNITY_CACHE_DIR']
        if not isinstance(cache_dir, str):
            cache_dir = cache_dir.abspath

        os.makedirs(cache_dir, exist_ok=True)
        target_base_name = os.path.basename(target)

        for idx in range(num_unity_files):
            unity_filename = f'{cache_dir}/{target_base_name}_{idx}.cpp'
            unity_source_files.append(unity_filename)
            begin = sources_per_file*idx
            end = sources_per_file*(idx+1)
            env.UnitySource(
                target = unity_filename,
                source = source_files[begin:end]
            )
        
        if len(other_nodes) > 0:
            print(f'Exluded {len(other_nodes)} node(s) from Unity build.')
        return [base_generator(target = target, source = unity_source_files + other_nodes, *args, **kwargs)]
    return generator

def _flatten_source(source : list):
    source_files = []
    other_nodes = []
    for ele in source:
        if isinstance(ele, list):
            more_sources, more_other = _flatten_source(ele)
            source_files.extend(more_sources)
            other_nodes.extend(more_other)
        elif isinstance(ele, File):
            source_files.append(ele.abspath)
        elif isinstance(ele, str):
            source_files.append(ele)
        else:
            other_nodes.append(ele)

    return source_files, other_nodes

def _generate_unity_file_msg(target, source, env : Environment):
    assert(len(target) == 1)
    return f'Generating {str(target[0])} from {len(source)} source files.'

def _generate_unity_file(target, source, env : Environment):
    assert(len(target) == 1)

    unity_filename = target[0].abspath
    with open(unity_filename, 'w') as f:
        for source_file in source:
            fpath = source_file.abspath.replace("\\", "\\\\")
            f.write(f'#include "{fpath}"\n')