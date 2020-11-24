#!/bin/bash

echo "Preparing the system to install..."
aclocal --force
libtoolize -c -f
automake -a -c -f
autoconf -i -f -v
echo "Done. Run the configure script to set up the makefiles"
