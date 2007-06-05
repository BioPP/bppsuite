#! /bin/sh
touch NEWS README AUTHORS ChangeLog
aclocal
libtoolize --copy
autoconf
automake --add-missing --copy
