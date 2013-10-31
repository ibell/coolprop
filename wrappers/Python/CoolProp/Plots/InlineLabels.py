# -*- coding: utf-8 -*-
"""
Created on Wed Oct 30 16:21:38 2013

@author: logan
"""

import math
import numpy
import matplotlib

from numpy.linalg import inv
from matplotlib import pyplot


def get_inline_label_bounds(theta, mat_text, axis=None):
    if axis is None:
        axis = pyplot.gca()
    renderer = axis.figure.canvas.get_renderer()

    # solve for the text width and height
    bbox = mat_text.get_window_extent(renderer)
    theta_mat = numpy.array([[math.cos(theta), math.sin(theta)],
                             [math.sin(theta), math.cos(theta)]])
    txt_size = numpy.dot(inv(theta_mat), numpy.array([bbox.width,
                                                      bbox.height]))

    # calculate the bounds (for the rotated bbox)
    bounds = bbox.corners()
    bounds[0][0] += txt_size[1] * math.sin(theta)
    bounds[1][1] -= txt_size[0] * math.sin(theta)
    bounds[2][1] += txt_size[0] * math.sin(theta)
    bounds[3][0] -= txt_size[1] * math.sin(theta)

    #print 'a'
    #print ind
    #print bbox.corners()
    #print bounds
    return bounds


def point_in_bbox(test_point, bounds):
    def distance(point1, point2):
        return numpy.sum((point2 - point1) ** 2) ** 0.5

    def perp(a):
        b = numpy.empty_like(a)
        b[0] = -a[1]
        b[1] = a[0]
        return b

    def seg_intersect(seg1, seg2):
        line_vector1 = seg1[1] - seg1[0]
        line_vector2 = seg2[1] - seg2[0]
        line_vector3 = seg1[0] - seg2[0]

        dap = perp(line_vector1)
        denom = numpy.dot(dap, line_vector2)
        numer = numpy.dot(dap, line_vector3)
        return (numer / denom) * line_vector2 + seg2[0]

    center_point = seg_intersect([bounds[0], bounds[3]],
                                 [bounds[1], bounds[2]])
    dist_center = numpy.array([distance(x, center_point) for x in bounds])
    dist_point =  numpy.array([distance(x, test_point) for x in bounds])
    return numpy.any((dist_point - dist_center) < 0.0)


class InlineLabel(object):
    def __init__(self, line, xval=None, yval=None, axis=None, **kwargs):
        self.line = line
        self.axis = axis
        if axis is None:
            self.axis = pyplot.gca()

        self.x_vals = self.line['x']
        self.y_vals = self.line['y']
        if len(self.x_vals) != len(self.y_vals):
            raise IndexError('asdasd')

        if xval is None and yval is None:
            raise ValueError("asdasdasd")

        if xval is not None:
            self.swap = False
            self.index = numpy.argmin(numpy.abs(self.x_vals - xval))
        elif yval is not None:
            self.swap = True
            self.index = numpy.argmin(numpy.abs(self.y_vals - yval))

        self.kwargs = kwargs
        self.kwargs.pop('verticalalignment', None)
        self.kwargs.pop('horizontalalignment', None)

        self.mat_text = None
        self.theta = None
        self.badness = []

    def calc_rotation(self, ind, radians=True):
        a = 1
        b = 0
        if ind == len(self.x_vals)-1:
            a = 0
            b = -1

        angle_trans = self.axis.transData.transform_angles
        theta = math.atan((self.y_vals[ind + a] - self.y_vals[ind + b]) /
                          (self.x_vals[ind + a] - self.x_vals[ind + b]))

        point = numpy.array([self.x_vals[ind], self.y_vals[ind]])
        theta = angle_trans(numpy.array((theta, )),
                            point.reshape((1, 2)),
                            radians=True)
        if radians:
            self.theta = float(theta)
        else:
            self.theta = float(math.degrees(theta))
        return self.theta

    def create_inline_label(self, ind=None):
        if ind is None:
            ind = 0
        xloc = self.x_vals[ind]
        yloc = self.y_vals[ind]
        self.calc_rotation(ind, radians=False)

        mat_text = pyplot.text(xloc, yloc,
                               self.line['label'],
                               verticalalignment='center',
                               horizontalalignment='center',
                               rotation=self.theta,
                               **self.kwargs)

        mat_text.set_rotation_mode('anchor')
        self.mat_text = mat_text
        return mat_text

    def adjust_inline_label(self, ind):
        self.calc_rotation(ind, radians=False)
        self.mat_text.set_rotation(self.theta)
        self.mat_text.set_position((self.x_vals[ind], self.y_vals[ind]))

        renderer = self.axis.figure.canvas.get_renderer()
        self.mat_text.draw(renderer)

    def get_overlap_badness(self, ind, labels=[]):
        if not labels:
            return 0

        badness = 0
        bounds = get_inline_label_bounds(self.theta, self.mat_text, self.axis)
        for label in labels:
            if isinstance(label, matplotlib.text.Text):
                rotation = label.get_rotation()
            elif isinstance(label, InlineLabel):
                rotation = label.theta
                label = label.mat_text

            other_bounds = get_inline_label_bounds(rotation, label, self.axis)
            for point in other_bounds:
                if point_in_bbox(point, bounds):
                    badness += 100
                    break

        return badness

    def get_satline_badness(self, ind, sat_lines=[]):
        if not sat_lines:
            return 0

        badness = 0
        bounds = get_inline_label_bounds(self.theta, self.mat_text, self.axis)
        for line in sat_lines:
            if line['type'] != 'Q':
                continue

            for (x, y) in zip(line['x'], line['y']):
                point = self.axis.transData.transform([x, y])
                if point_in_bbox(point, bounds):
                    badness += 100
                    break

        return badness

    def get_axis_limit_badness(self, ind):
        axis_trans = self.axis.transData
        xmin, xmax = axis_trans.transform(self.axis.get_xlim())
        ymin, ymax = axis_trans.transform(self.axis.get_ylim())

        badness = 0
        bounds = get_inline_label_bounds(self.theta, self.mat_text, self.axis)
        for point in bounds:
            #print ind, point, [xmin, ymin, xmax, ymax]
            if point[0] < xmin or point[0] > xmax:
                badness += 100
                break
            elif point[1] < ymin or point[1] > ymax:
                badness += 100
                break
        return badness

    def get_badness(self, sat_lines=[], labels=[]):
        self.badness = [abs(i - self.index) for i in range(len(self.x_vals))]
        self.badness = numpy.array(self.badness, dtype=float)

        for i, val in enumerate(self.x_vals):
            self.adjust_inline_label(i)
            self.badness[i] += self.get_satline_badness(i, sat_lines)
            #self.badness[i] += self.get_overlap_badness(i, labels)
            #self.badness[i] += self.get_axis_limit_badness(i)
        #print
        #print self.badness

    def place_label(self, method='auto', labels=[], sat_lines=[]):
        axis_autoscale = self.axis.get_autoscalex_on()
        self.axis.set_autoscalex_on(False)

        if method == 'center' or method == 'manual':
            ind = numpy.argmin(self.badness)
            self.create_inline_label(ind)
            self.index = ind
        elif method == 'auto':
            self.create_inline_label()
            self.get_badness(sat_lines, labels)
            ind = numpy.argmin(self.badness)
            self.adjust_inline_label(ind)
            self.index = ind
        else:
            raise KeyError('asdasd')

        if axis_autoscale:
            self.axis.set_autoscalex_on(True)
        return self.mat_text


def main():
    line = {'x': numpy.array([1, 2, 3, 4, 5]),
            'y': numpy.array([1, 2, 3, 4, 5]),
            'label': '12'}

    pyplot.plot(line['x'], line['y'])
    text = pyplot.text(3, 3, 'FF')
    #pyplot.axis([0, 7, 0, 7])

    label = InlineLabel(line, xval=4)
    label.place_label(labels=[text])

    pyplot.show()


if __name__ == "__main__":
    main()