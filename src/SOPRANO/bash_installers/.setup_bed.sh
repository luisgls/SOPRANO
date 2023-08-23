#!/bin/bash

# Build bedtools
if [ ! -f "$BEDTOOLS_EXEC_PATH" ]
then
	cd "$TOOLS_DIR_PATH" || exit
	echo "-- downloading bedtools binary..."
	wget https://github.com/arq5x/bedtools2/releases/download/v2.31.0/bedtools.static
	mv bedtools.static "$BEDTOOLS_EXEC_PATH"
	chmod a+x "$BEDTOOLS_EXEC_PATH"
	echo "-- bedtools executable configured"
	cd "$PROJECT_ROOT" || exit
else
	echo "-- bedtools detected"
fi