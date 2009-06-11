#! /bin/sh
touch NEWS README AUTHORS ChangeLog
aclocal
libtoolize --copy
autoconf
automake --add-missing --copy
#Hack for the AC_CONFIG_MACRO_DIR([m4]):
if [ ! -d m4 ]; then mkdir m4; fi
