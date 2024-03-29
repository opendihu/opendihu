"""SCons.Tool.ilink32

XXX

"""

#
# __COPYRIGHT__
#
# Permission is hereby granted, free of charge, to any person obtaining
# a copy of this software and associated documentation files (the
# "Software"), to deal in the Software without restriction, including
# without limitation the rights to use, copy, modify, merge, publish,
# distribute, sublicense, and/or sell copies of the Software, and to
# permit persons to whom the Software is furnished to do so, subject to
# the following conditions:
#
# The above copyright notice and this permission notice shall be included
# in all copies or substantial portions of the Software.
#
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY
# KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE
# WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
# NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE
# LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION
# OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION
# WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
#

__revision__ = "__FILE__ __REVISION__ __DATE__ __DEVELOPER__"

import SCons.Tool
import SCons.Tool.bcc32
import SCons.Util

def generate(env) -> None:
    """Add Builders and construction variables for Borland ilink to an
    Environment."""
    SCons.Tool.createSharedLibBuilder(env)
    SCons.Tool.createProgBuilder(env)

    env['LINK']        = '$CC'
    env['LINKFLAGS']   = SCons.Util.CLVar('')
    env['LINKCOM']     = '$LINK -q $LINKFLAGS -e$TARGET $SOURCES $LIBS'
    env['LIBDIRPREFIX']=''
    env['LIBDIRSUFFIX']=''
    env['LIBLINKPREFIX']=''
    env['LIBLINKSUFFIX']='$LIBSUFFIX'


def exists(env):
    # Uses bcc32 to do linking as it generally knows where the standard
    # LIBS are and set up the linking correctly
    return SCons.Tool.bcc32.findIt('bcc32', env)

# Local Variables:
# tab-width:4
# indent-tabs-mode:nil
# End:
# vim: set expandtab tabstop=4 shiftwidth=4:
