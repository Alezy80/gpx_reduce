#!/usr/bin/env python
# -*- coding: utf8 -*-

'''
gpx_reduce: removes points from gpx-files to reduce filesize and
tries to keep introduced distortions to the track at a minimum.
Copyright (C) 2011 travelling_salesman on OpenStreetMap
 
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
'''


import sys
import pylab as pl
import scipy as sc
from math import *
from lxml import etree
from optparse import OptionParser


parser = OptionParser('usage: %prog [options] input-file.gpx')
parser.add_option('-p', '--plot', action='store_true', dest='plot',
    default=False, help='Show a plot of the result at the end.')
parser.add_option('-d', '--dist', action='store', type='float', dest='max_dist',
    default=1.0, help='Maximum distance of line from original points in meters')
parser.add_option('-o', '--out', action='store', type='string',
    dest='ofname', default=None, help='Output file name')
parser.add_option('-m', '--maxsep', action='store', type='float', dest='max_sep',
    default=200.0, help='Maximum separation of points. No points will be deleted where the resulting distance would become greater than maxsep. Standard JOSM settings will not display points spaced more than 200m. Zero value means no limit.')
(options, args) = parser.parse_args()

if len(args) < 1:
    parser.print_usage()
    exit(2)


# use the WGS-84 ellipsoid
rE = 6356752.314245 # earth's radius
a = 6378137.0
b = 6356752.314245179


def norm(v):
    return sqrt(sum([i**2 for i in v]))

def linedistance_position(p1, pm, p2):
    # returns distance of pm from line between p1 and p2
    line = p2 - p1
    linel = norm(line)
    vm = pm - p1
    if linel == 0.0:
        return norm(vm), 0.5
    linem = line / linel
    
    position = pl.dot(vm, linem) / linel
    distance = vm - line * position
    
    return norm(distance), position

def get_xyz(dicti):
    return dicti['x'], dicti['y'], dicti['z']

def rotate(x, y, phi):
    return x*cos(phi) - y*sin(phi), x*sin(phi) + y*cos(phi)
    
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
    r = ele + a * b / norm([s*a, c*b])
    lon = radians(lon)
    return r * c * sin(lon), r * c * (-cos(lon)), r * s

def xyz_to_latlonele(x, y, z):
    r = norm([x, y, z])
    if (r == 0):
        return 0.0, 0.0, 0.0
    lat = degrees(atan2(z, norm([x, y])))
    lon = degrees(atan2(x, -y))
    ele = r * (1.0 - a * b / norm([a*z, b*x, b*y]))
    return lat, lon, ele



for fname in args:
    # initialisations
    tracksegs_old = []
    tracksegs_new = []
    sumx = 0.0
    sumy = 0.0
    sumz = 0.0
    
    # import xml data from files
    print 'opening file', fname
    infile = open(fname)
    
    tree = etree.parse(infile)
    infile.close()
    gpx = tree.getroot()
    if gpx.nsmap.has_key(None):
        nsmap = '{' + gpx.nsmap[None] + '}'
    else:
        nsmap = ''
                
    
    # extract data from xml
    for trkseg in gpx.findall('.//' + nsmap + 'trkseg'):
        trkpts = trkseg.findall(nsmap + 'trkpt')
        n = len(trkpts)
        
        # extract coordinate values
        lats = [float(trkpt.get('lat')) for trkpt in trkpts]
        lons = [float(trkpt.get('lon')) for trkpt in trkpts]
        eles = [float(trkpt.find(nsmap + 'ele').text) for trkpt in trkpts]
        
        # save original trackseg for plotting
        if options.plot:
            tracksegs_old.append([[lats[i], lons[i], eles[i]] for i in range(n)])
        
        # calculate projected points to work on
        points = [{} for i in range(n)]
        for i in range(n):
            x, y, z = latlonele_to_xyz(lats[i], lons[i], eles[i])
            points[i]['x'] = x
            points[i]['y'] = y
            points[i]['z'] = z
            sumx += x
            sumy += y
            sumz += z
        
        # create lists of connections to all previous points
        # and distances to intermediate points
        points[0]['distances'] = {}
        for i2 in range(1, n):
            points[i2]['distances'] = {i2-1:0.0}
            for i1 in reversed(range(i2-1)):
                p1 = sc.array(get_xyz(points[i1]))
                p2 = sc.array(get_xyz(points[i2]))
                if 0.0 < options.max_sep and options.max_sep <= norm(p2 - p1):
                    break # point separation is too far
                
                ok = True
                dlist = []
                # go through range(i1+1, i2) but start in the middle
                for im in range((i1-1+i2)/2, i1, -1) + range((i1+1+i2)/2, i2):
                    pm = sc.array(get_xyz(points[im]))
                    d, l = linedistance_position(p1, pm, p2)
                    if (l >= 0.0 and l <= 1.0 and d <= options.max_dist):
                        dlist.append(d)
                    else:
                        ok = False
                        break
                if ok:
                    points[i2]['distances'][i1] = sum(
                        [(i / options.max_dist)**2 for i in dlist])
        
        # execute routing algorithm on points
        points[0]['cost'] = 1.0
        points[0]['prev'] = -1
        for i in range(1, n):
            imin = None
            costmin = float('inf')
            for prev, dist in (points[i]['distances']).iteritems():
                cost = points[prev]['cost'] + 1.0 + dist
                if cost < costmin:
                    imin = prev
                    costmin = cost
            points[i]['cost'] = costmin
            points[i]['prev'] = imin
        
        # trace route backwards to collect final points
        final_pnums = []
        i = n-1
        while i >= 0:
            final_pnums = [i] + final_pnums
            i = points[i]['prev']
        
        n_new = len(final_pnums)
        print 'number of points:', n, '-', n-n_new, '=', n_new
        
        # delete certain points from original data
        delete_pnums = [i for i in range(n) if i not in final_pnums]
        for i in reversed(delete_pnums):
            del trkseg[trkseg.index(trkpts[i])]
        
        # save reduced trackseg for plotting
        if options.plot:
            tracksegs_new.append([
                [float(trkpt.get('lat')), float(trkpt.get('lon')), float(trkpt.find(nsmap + 'ele').text)]
                for trkpt in trkseg.findall(nsmap + 'trkpt')])
    
        
    # export data to file
    if options.ofname != None:
        ofname = options.ofname
    elif fname.endswith('.gpx'):
        ofname = fname[:-4] + '_reduced.gpx'
    else:
        ofname = fname + '_reduced.gpx'
    outfile = open(ofname, 'w')
    outfile.write(etree.tostring(tree, xml_declaration=True,
        pretty_print=True, encoding='utf-8'))
    outfile.close()
    print 'modified copy written to', ofname
    
    
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