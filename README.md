# gpx reduce

For more information see the [original script](http://wiki.openstreetmap.org/wiki/User:Travelling_salesman/gpx_reduce).

This is a python script that removes (almost) unnecessary points from gpx-files. Set your gps receiver to the highest point density available (e.g. 1s) and use this script on your gpx-files before osm-upload. The script yields highly optimized results.

The script must be saved to a file (gpx_reduce.py) and the executable flag should be set. Then it can be executed from the command line. For more documentation, type 
```sh
./gpx_reduce.py --help
```

If the script produces an output that looks non-optimal, this might be because of the third dimension. Otherwise leave me a message.

### Requirements
* [python](http://en.wikipedia.org/wiki/Python_(programming_language)) interpreter
* [scipy](http://www.scipy.org/) scientific tools for python
* [lxml](http://lxml.de/) xml library
* matplotlib

On a debian based linux distribution (such as Ubuntu) you can install these with:
```sh
sudo apt-get install python python-scipy python-lxml python-matplotlib 
```

### Algorithm

gpx_reduce finds the track that globally optimizes a given cost function under the restriction that no original point may be further away from the resulting track than a certain limit. Points can be deleted but will never be moved or inserted.

The cost function can be easily modified. In the standard setting it consists of the number of remaining points plus (normalized) squared distances to removed points, plus a penalty for large-angle kinks at the remaining trackpoints.

The number of possible tracks is huge and not all of them could be explored individually to find the optimum. Instead the a modified version of the Dijkstra algorithm is utilized. It successively finds the shortest route to each point and then only has to trace back to all possible predecessor points instead of routes. Therefore the best trace is found and the task is performed in a reasonable amount of time. However if the track is very long and offers many possibilities, it may still take several seconds to compute.

Points will never become separated beyond a certain limit (200m by default), since josm would not show lines between further separated points with default settings. There is an additional speed-dependent distance-limit that keeps points more evenly spaced. 

### Code

See in repo

### License
```
This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program.  If not, see <http://www.gnu.org/licenses/>.
```