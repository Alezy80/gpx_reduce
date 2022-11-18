# gpx reduce

For more information see the [original script](http://wiki.openstreetmap.org/wiki/User:Travelling_salesman/gpx_reduce).

This is a python script that removes (almost) unnecessary points from gpx-files. Set your gps receiver to the highest point density available (e.g. 1s) and use this script on your gpx-files before osm-upload. The script yields highly optimized results.

The script must be saved to a file (gpx_reduce.py) and the executable flag should be set. Then it can be executed from the command line. For more documentation, type 
```sh
./gpx_reduce.py --help
```

If the script produces an output that looks non-optimal, this might be because of the third dimension. Otherwise leave me a message.

### Changes from original version
For now I don't changed any algorithms in the script, but only added some features:
* finetuned output abilities, like adding more compact output.
* allowing to strip some unwanted tag (i.e. you can anonymize your track time and speed)
* more robust handling files without elevation
* time parsing conforms to iso8601
* not only tracks can be reduced, but also routes
* supports globs in file name
* I've made some changes towards make porting to Python 3 more easily

### Requirements
* [python](http://en.wikipedia.org/wiki/Python_(programming_language)) interpreter
* [numpy](https://numpy.org/) the fundamental package for scientific computing with Python
* [matplotlib](http://matplotlib.org/) 2D plotting library
* [lxml](https://pypi.python.org/pypi/lxml) xml library
* [iso8601](http://pypi.python.org/pypi/iso8601) date/time library

On a Debian based Linux distribution (such as Ubuntu) you can install these with:
```sh
sudo apt-get install python python-numpy python-lxml python-matplotlib python-iso8601
```

For Windows users, copy the requirements.txt and run this in your `cmd` prompt :
```sh
pip install -r requirements.txt
```

### Algorithm

gpx_reduce finds the track that globally optimizes a given cost function under the restriction that no original point may be further away from the resulting track than a certain limit. Points can be deleted but will never be moved or inserted.

The cost function can be easily modified. In the standard setting it consists of the number of remaining points plus (normalized) squared distances to removed points, plus a penalty for large-angle kinks at the remaining trackpoints.

The number of possible tracks is huge and not all of them could be explored individually to find the optimum. Instead the a modified version of the Dijkstra algorithm is utilized. It successively finds the shortest route to each point and then only has to trace back to all possible predecessor points instead of routes. Therefore the best trace is found and the task is performed in a reasonable amount of time. However if the track is very long and offers many possibilities, it may still take several seconds to compute.

Points will never become separated beyond a certain limit (200m by default), since josm would not show lines between further separated points with default settings. There is an additional speed-dependent distance-limit that keeps points more evenly spaced. 

### Code

See in repo

### Example

```sh
./gpx_reduce.py -d 10 -n 60 -t 60 -m 1000 -p -2 my_track.gpx
```
Removes unnecessary points making new track not farther that 10 meters from original. It you stay still the points are written not more often than once in 60 seconds. The track segments can be up to 1000 meters. The only Latitude and Longitude is retained. In my opinion theese settings are good compromise between accuracy and resulting file size for tracks got by driving, but if you want to draw some footways in forest, it's better to use more conservative settings.

For a bike track with an Ultratrac gps mode (Garmin), this seems to work well:
```sh
./gpx_reduce.py -d 600 -m 600 -t 50 -4 my_track.gpx
```

where `-4` keeps the elevation values (check help).


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