#!/bin/sh

package_dir="osu/megraw/llscan"
package="osu.megraw.llscan"
main="Scan"

cliCP="../lib/commons-cli-1.3.1.jar"

ver=`javac -version`
echo $ver

javac -cp $cliCP $package_dir/*.java
jar cvfe scan.jar $package.$main $package_dir/*.class
rm $package_dir/*.class
mv scan.jar ../main
