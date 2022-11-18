#!/usr/bin/env python2

"""
gpx_reduce v1.8: removes points from gpx-files to reduce filesize and
tries to keep introduced distortions to the track at a minimum.
Copyright (C) 2011,2012,2013,2015,2016,2017 travelling_salesman on OpenStreetMap

changelog: v1.2: clarity refractoring + speedup for identical points
           v1.3: new track weighting functions, progress display
           v1.4: new track weighting function, restructuring for memory saving
           v1.5: algorithm speedup by roughly a factor of 2 by eliminating some cases.
           v1.6: presets for train etc.
           v1.7: introduced weighting function for elevation errors
           v1.8: speed-dependent distance limit

Copyright (C) 2017,2018 Alezy at github.com/Alezy80
changelog: v1.8.1: adding compact output, allows to strip unneeded tags
           v1.8.2: supports processing of tracks and routes, globs in filenames
                   under the Windows too, more robust working, tuned strip
                   options
           v1.8.3: Faster processing (about 35%).

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
"""

from __future__ import print_function

import glob
import time
from math import *
from optparse import OptionParser, OptionGroup
from sys import stdout

import pylab as pl
from iso8601 import parse_date
from lxml import etree
from numpy import array, dot

parser = OptionParser('usage: %prog [options] input-file.gpx')
parser.add_option('-v', '--verbose', action='store', type='int',
                  dest='verbose', default=1, help='verbose=[0,1]')
parser.add_option('-p', '--plot', action='store_true', dest='plot',
                  default=False, help='Show a plot of the result at the end.')
parser.add_option('-d', '--dist', action='store', type='float', dest='max_dist',
                  default=1.0, help='Maximum distance of line from original points in meters')
parser.add_option('-o', '--out', action='store', type='string',
                  dest='ofname', default=None, help='Output file name')
parser.add_option('-m', '--maxsep', action='store', type='float', dest='max_sep',
                  default=200.0,
                  help='Absolute maximum separation of points. No points will be deleted where the resulting distance '
                       'would become greater than maxsep. Standard JOSM settings will not display points spaced more '
                       'than 200m. -1 means no limit.')
parser.add_option('-n', '--maxsep0', action='store', type='float', dest='max_sep0',
                  default=10.0, help='Maximum separation of points at zero speed.')
parser.add_option('-t', '--maxseptime', action='store', type='float', dest='max_sep_time',
                  default=10.0,
                  help='Maximum time separation of points, which will be added to maxsep0. Maximum allowed point '
                       'separation will be min(m, n + t*v) where v is the speed in m/s.')
parser.add_option('-e', '--ele-weight', action='store', type='float',
                  dest='ele_weight', default=0.0,
                  help='Weighting of elevation errors vs. horizontal errors. Default is 0.')
parser.add_option('-b', '--bend', action='store', type='float', dest='bend',
                  default=1.0,
                  help='Penalty value for large bending angles at each trackpoint. Larger values (1 or more) make the '
                       'track smoother.')
parser.add_option('-w', '--weighting', action='store', type='string',
                  dest='weighting', default='exp',
                  help='''Weighting function to be minimized for track reduction:
pnum (number of points),
sqrdistsum (number of points plus sum of squared distances to leftout points),
sqrdistmax (number of points plus sum of squared distances to each maximally separated leftout point per
new line segment),
sqrlength (number of points plus sum of squared new line segment lengths normalized by maxsep),
mix (number of points plus sum of squared distances to each maximally separated leftout point per new line segment
weighted with corresponding segment length),
exp (number of points plus sum of squared distances to leftout points with exponential weighting of 1/2, 1/4, 1/8...
from furthest to closest point). exp=standard''')
strip_group = OptionGroup(parser,
                          'Strip sensitive info, '
                          'You can remove unwanted info from trackpoints, but it does not remove any metadata! All '
                          'options will always retain latitude and longitude')
strip_group.add_option('-c', '--compact', action='store_true', dest='compact_output',
                       default=False, help='Makes output file more compact by removing spaces (pretty_print=false)')
strip_group.add_option('-2', action='store_const', dest='keep_tags', const='',
                       help='Retain only latitude, longitude')
strip_group.add_option('-3', action='store_const', dest='keep_tags', const='ele',
                       help='Retain only latitude, longitude, elevation')
strip_group.add_option('-4', action='store_const', dest='keep_tags', const='ele,time',
                       help='Retain only latitude, longitude, ele and time')
strip_group.add_option('-s', '--strip', action='store', type='string', dest='keep_tags',
                       default='***',
                       help='Strip all track point informaton EXCEPT listed here tags, splitted by comma. Generally, '
                            'you may want to keep elevation (ele) or time (time) tags. If you want to keep it both, '
                            'write option in this way "-s ele,time"')
parser.add_option_group(strip_group)

(options, args) = parser.parse_args()

if len(args) < 1:
    parser.print_usage()
    exit(2)

# use the WGS-84 ellipsoid
rE = 6356752.314245  # earth's radius
a = 6378137.0
b = 6356752.314245179
reduced_suffix = '_reduced.gpx'


def norm(coords):
    if len(coords) == 2:
        return sqrt(coords[0] * coords[0] + coords[1] * coords[1])
    if len(coords) == 3:
        return sqrt(coords[0] * coords[0] + coords[1] * coords[1] + coords[2] * coords[2])

    result = 0  # not pythonic way, but much faster than iterator
    for coord in coords:
        result += coord * coord
    return sqrt(result)


def distance(p1_, pm_, p2_, ele_weight=1.0):
    # returns distance of pm from line between p1 and p2

    p1, pm, p2 = array(p1_), array(pm_), array(p2_)
    h1, hm, h2 = norm(p1), norm(pm), norm(p2)
    if ele_weight != 1.0 and min(h1, hm, h2) > 0.0:
        hmean = (h1 + hm + h2) / 3.0
        p1 *= (ele_weight + (1.0 - ele_weight) * hmean / h1)
        pm *= (ele_weight + (1.0 - ele_weight) * hmean / hm)
        p2 *= (ele_weight + (1.0 - ele_weight) * hmean / h2)
    line = p2 - p1
    linel = norm(line)
    vm = pm - p1
    if linel == 0.0:
        return norm(vm)
    linem = line / linel

    position = pl.dot(vm, linem) / linel
    if position < 0.0:
        return norm(vm)
    elif position > 1.0:
        return norm(pm - p2)
    else:
        return norm(vm - line * position)


def rotate(x, y, phi):
    return x * cos(phi) - y * sin(phi), x * sin(phi) + y * cos(phi)


def project_to_meters(lat, lon, latm, lonm):
    # azimuthal map projection centered at average track coordinate
    lon -= lonm
    xyz = latlonele_to_xyz(lat, lon, 0.0)
    zy = rotate(xyz[2], xyz[1], radians(90 - latm))
    lat2 = atan2(zy[0], norm([zy[1], xyz[0]]))
    lon2 = atan2(xyz[0], -zy[1])
    x_meters = rE * sin(lon2) * (pi / 2.0 - lat2)
    y_meters = -rE * cos(lon2) * (pi / 2.0 - lat2)
    return x_meters, y_meters


def latlonele_to_xyz(lat, lon, ele):
    s = sin(radians(lat))
    c = cos(radians(lat))
    r = ele + a * b / norm([s * a, c * b])
    lon = radians(lon)
    return r * c * sin(lon), r * c * (-cos(lon)), r * s


def xyz_to_latlonele(x, y, z):
    r = norm([x, y, z])
    if r == 0:
        return 0.0, 0.0, 0.0
    lat = degrees(atan2(z, norm([x, y])))
    lon = degrees(atan2(x, -y))
    ele = r * (1.0 - a * b / norm([a * z, b * x, b * y]))
    return lat, lon, ele


def reduced_track_indices(coordinate_list, timesteps=None):
    # returns a list of indices of trackpoints that constitute the reduced track
    # takes a list of kartesian coordinate tuples
    m = len(coordinate_list)
    if m == 0:
        return []
    if timesteps is not None and len(timesteps) != len(coordinate_list):
        timesteps = None

    # number of dimensions
    d = len(coordinate_list[0])

    # remove identical entries (can speed up algorithm considerably)
    original_indices = [0]
    points = [{'p': coordinate_list[0], 'weight': 1}]
    if timesteps is not None:
        points[0]['t'] = timesteps[0]
    for i in range(1, m):
        if False in [coordinate_list[i - 1][j] == coordinate_list[i][j] for j in range(d)]:
            original_indices.append(i)
            points.append({'p': coordinate_list[i], 'weight': 1})
            if timesteps is not None:
                points[-1]['t'] = timesteps[i]
        else:
            points[-1]['weight'] += 1
    n = len(points)

    # progress printing initialisations
    progress_printed = False
    progress = None
    tprint = time.time()

    # execute Dijkstra-like algorithm on points
    points[0]['cost'] = 1.0
    points[0]['prev'] = -1

    for i2 in range(1, n):
        penalties = {}
        costmin = float('inf')
        for i1 in reversed(list(range(i2))):
            p1 = array(points[i1]['p'])
            p2 = array(points[i2]['p'])
            seglength = norm(p2 - p1)

            # estimate speed between p1 and p2
            if timesteps is not None:
                dt = (points[i2]['t'] - points[i1]['t']).total_seconds()
                v = seglength / max(0.1, dt)
            else:
                v = seglength / float(i2 - i1)  # assume 1s time spacing

            max_sep = options.max_sep0 + v * options.max_sep_time
            if options.max_dist >= 0:
                max_sep = min(max_sep, options.max_sep)

            if seglength >= max_sep and i1 != i2 - 1:
                # point separation is too far
                # but always accept direct predecessor i1 = i2 - 1
                if seglength >= max_sep + options.max_dist:
                    # no chance to find a valid earlier predecessor point
                    break
                else:
                    continue

            if points[i1]['cost'] + 1.0 > costmin:
                # the possible predecessor i1 is already too bad.
                continue

            i1_i2_segment_valid = True
            lower_i1_possible = True
            distance_squaremax = 0.0
            distance_squaresum = 0.0
            distances_squared = []
            # iterate all medium points between i1 and i2
            for im in range(i1 + 1, i2):
                pm = array(points[im]['p'])
                d = distance(p1, pm, p2, options.ele_weight)
                if d <= options.max_dist:
                    d_sq = (d / options.max_dist) ** 2
                    distance_squaremax = max(distance_squaremax, d_sq)
                    distance_squaresum += points[im]['weight'] * d_sq
                    distances_squared.append(d_sq)
                else:
                    i1_i2_segment_valid = False

                    # check if connection to any further point i1 is impossible
                    d1 = pl.dot(p1 - p2, p1 - p2)
                    d2 = pl.dot(pm - p2, pm - p2)
                    dd = options.max_dist ** 2
                    d1d2 = pl.dot(p1 - p2, pm - p2)
                    # formula from cosines of point separation angle and cone-opening angles around points
                    if d1 > dd and d2 > dd and (d1d2 + dd) ** 2 < (d2 - dd) * (d1 - dd):
                        lower_i1_possible = False
                        break

            if not lower_i1_possible:
                break

            if i1_i2_segment_valid:
                if options.weighting == 'sqrdistmax':
                    penalties[i1] = distance_squaremax
                elif options.weighting == 'sqrdistsum':
                    penalties[i1] = distance_squaresum
                elif options.weighting == 'sqrlength':
                    penalties[i1] = (seglength / max_sep) ** 2
                elif options.weighting == 'mix':
                    penalties[i1] = (distance_squaremax * (1.0 + seglength / max_sep))
                elif options.weighting == 'exp':
                    penalties[i1] = 0.5 * sum([0.5 ** i * d for i, d in
                                               enumerate(sorted(distances_squared, reverse=True))])
                else:
                    penalties[i1] = 0.0

                # add a penalty for kinks
                if options.bend > 0.:
                    if points[i1]['prev'] != -1:
                        p0 = array(points[points[i1]['prev']]['p'])
                        v0 = p1 - p0
                        v1 = p2 - p1
                        if norm(v0) > 0. and norm(v1) > 0.:
                            v0 /= norm(v0)
                            v1 /= norm(v1)
                            kink = (1.0 - dot(v0, v1)) / 2.0
                            penalties[i1] += options.bend * kink

        # find best predecessor
        imin = None
        costmin = float('inf')
        for prev, penalty in penalties.items():
            # cost function is sum of points used (1.0) plus penalties
            cost = points[prev]['cost'] + 1.0 + penalty
            if cost < costmin:
                imin = prev
                costmin = cost
        points[i2]['cost'] = costmin
        points[i2]['prev'] = imin

        # print progess
        # if options.verbose == 1 and (100 * i2) / n > progress and time.time() >= tprint + 1:
            # tprint = time.time()
            # progress = (100 * i2) / n
            # print('\r', progress, '% of', n, 'points', end='')
            # stdout.flush()
            # progress_printed = True

    if progress_printed:
        print('\r', end='')

    # trace route backwards to collect final points
    final_pnums = []
    i = n - 1
    while i >= 0:
        final_pnums = [i] + final_pnums
        i = points[i]['prev']

    return [original_indices[i] for i in final_pnums]


def strip_trackpoint_tags(trkpts, nsmap):
    if options.keep_tags == '***':
        return

    keep_tag_names = [nsmap + tag_name.strip() for tag_name in options.keep_tags.split(',')]

    for trkpt in trkpts:
        for child in trkpt:
            if child.tag not in keep_tag_names:
                trkpt.remove(child)


def get_ele(trkpt, nsmap):
    try:
        return float(trkpt.find(nsmap + 'ele').text)
    except AttributeError:
        return 0.0


def process_file(fname):
    # initialisations
    tracksegs_old = []
    tracksegs_new = []
    sumx, sumy, sumz = 0.0, 0.0, 0.0

    # import xml data from files
    print('opening file', fname)
    infile = open(fname)

    eparser = etree.XMLParser(remove_blank_text=True)
    tree = etree.parse(infile, eparser)
    infile.close()
    gpx = tree.getroot()
    if None in gpx.nsmap:
        nsmap = '{' + gpx.nsmap[None] + '}'
    else:
        nsmap = ''

    # extract data from xml. Process both tracks and routes
    for seg_type, pt_type in (('trkseg', 'trkpt'), ('rte', 'rtept')):
        for trkseg in gpx.findall('.//' + nsmap + seg_type):
            trkpts = trkseg.findall(nsmap + pt_type)
            n = len(trkpts)

            # extract coordinate values
            lats = [float(trkpt.get('lat')) for trkpt in trkpts]
            lons = [float(trkpt.get('lon')) for trkpt in trkpts]
            eles = [get_ele(trkpt, nsmap) for trkpt in trkpts]
            try:
                times = [parse_date(trkpt.find(nsmap + 'time').text) for trkpt in trkpts]
            except Exception as ex:
                print('There is no time in file.', ex)
                times = None

            # save original trackseg for plotting
            if options.plot:
                tracksegs_old.append([[lats[i], lons[i], eles[i]] for i in range(n)])

            # calculate projected points to work on
            coords = []
            for i in range(n):
                x, y, z = latlonele_to_xyz(lats[i], lons[i], eles[i])
                coords.append((x, y, z))
                sumx += x
                sumy += y
                sumz += z

            # execute the reduction algorithm
            final_pnums = reduced_track_indices(coords, times)

            n_new = len(final_pnums)
            print('number of points:', n, '-', n - n_new, '=', n_new)

            # delete certain points from original data
            delete_pnums = [i for i in range(n) if i not in final_pnums]
            for i in reversed(delete_pnums):
                del trkseg[trkseg.index(trkpts[i])]

            # save reduced trackseg for plotting
            if options.plot:
                tracksegs_new.append([[float(trkpt.get('lat')), float(trkpt.get('lon')), get_ele(trkpt, nsmap)]
                                      for trkpt in trkseg.findall(nsmap + pt_type)])

            strip_trackpoint_tags(trkpts, nsmap)

    # export data to file
    if options.ofname is not None:
        ofname = options.ofname
    elif fname.lower().endswith('.gpx'):
        ofname = fname[:-4] + reduced_suffix
    else:
        ofname = fname + reduced_suffix
    outfile = open(ofname, 'wb')
    outfile.write(etree.tostring(tree, xml_declaration=True,
                                 pretty_print=not options.compact_output, encoding='utf-8'))
    outfile.close()
    print('modified copy written to', ofname)

    # plot result to screen
    if options.plot:
        latm, lonm, elesum = xyz_to_latlonele(sumx, sumy, sumz)

        for trkseg in tracksegs_old:
            y_old = []
            x_old = []
            for trkpt in trkseg:
                xy = project_to_meters(trkpt[0], trkpt[1], latm, lonm)
                x_old.append(xy[0])
                y_old.append(xy[1])
            pl.plot(x_old, y_old, 'r.-')

        for trkseg in tracksegs_new:
            y_new = []
            x_new = []
            for trkpt in trkseg:
                xy = project_to_meters(trkpt[0], trkpt[1], latm, lonm)
                x_new.append(xy[0])
                y_new.append(xy[1])
            pl.plot(x_new, y_new, 'b.-')
        pl.grid()
        pl.gca().set_aspect('equal')
        pl.xlabel('x [m]')
        pl.ylabel('y [m]')
        pl.show()


for file_name in args:
    start_time = time.time()

    for file_name2 in glob.iglob(file_name):
        if reduced_suffix in file_name2 and '*' in file_name:
            continue
        process_file(file_name2)

    end_time = time.time()
    print('elapsed time=', end_time - start_time)
