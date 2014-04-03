(function() {
    !function() {
        var d3 = {
            version: "3.4.4"
        };
        var abs = Math.abs;
        d3.geo = {};
        d3.geo.stream = function(object, listener) {
            if (object && d3_geo_streamObjectType.hasOwnProperty(object.type)) {
                d3_geo_streamObjectType[object.type](object, listener);
            } else {
                d3_geo_streamGeometry(object, listener);
            }
        };
        function d3_geo_streamGeometry(geometry, listener) {
            if (geometry && d3_geo_streamGeometryType.hasOwnProperty(geometry.type)) {
                d3_geo_streamGeometryType[geometry.type](geometry, listener);
            }
        }
        var d3_geo_streamObjectType = {
            Feature: function(feature, listener) {
                d3_geo_streamGeometry(feature.geometry, listener);
            },
            FeatureCollection: function(object, listener) {
                var features = object.features, i = -1, n = features.length;
                while (++i < n) d3_geo_streamGeometry(features[i].geometry, listener);
            }
        };
        var d3_geo_streamGeometryType = {
            Sphere: function(object, listener) {
                listener.sphere();
            },
            Point: function(object, listener) {
                object = object.coordinates;
                listener.point(object[0], object[1], object[2]);
            },
            MultiPoint: function(object, listener) {
                var coordinates = object.coordinates, i = -1, n = coordinates.length;
                while (++i < n) object = coordinates[i], listener.point(object[0], object[1], object[2]);
            },
            LineString: function(object, listener) {
                d3_geo_streamLine(object.coordinates, listener, 0);
            },
            MultiLineString: function(object, listener) {
                var coordinates = object.coordinates, i = -1, n = coordinates.length;
                while (++i < n) d3_geo_streamLine(coordinates[i], listener, 0);
            },
            Polygon: function(object, listener) {
                d3_geo_streamPolygon(object.coordinates, listener);
            },
            MultiPolygon: function(object, listener) {
                var coordinates = object.coordinates, i = -1, n = coordinates.length;
                while (++i < n) d3_geo_streamPolygon(coordinates[i], listener);
            },
            GeometryCollection: function(object, listener) {
                var geometries = object.geometries, i = -1, n = geometries.length;
                while (++i < n) d3_geo_streamGeometry(geometries[i], listener);
            }
        };
        function d3_geo_streamLine(coordinates, listener, closed) {
            var i = -1, n = coordinates.length - closed, coordinate;
            listener.lineStart();
            while (++i < n) coordinate = coordinates[i], listener.point(coordinate[0], coordinate[1], coordinate[2]);
            listener.lineEnd();
        }
        function d3_geo_streamPolygon(coordinates, listener) {
            var i = -1, n = coordinates.length;
            listener.polygonStart();
            while (++i < n) d3_geo_streamLine(coordinates[i], listener, 1);
            listener.polygonEnd();
        }
        function d3_noop() {}
        function d3_adder() {}
        d3_adder.prototype = {
            s: 0,
            t: 0,
            add: function(y) {
                d3_adderSum(y, this.t, d3_adderTemp);
                d3_adderSum(d3_adderTemp.s, this.s, this);
                if (this.s) this.t += d3_adderTemp.t; else this.s = d3_adderTemp.t;
            },
            reset: function() {
                this.s = this.t = 0;
            },
            valueOf: function() {
                return this.s;
            }
        };
        var d3_adderTemp = new d3_adder();
        function d3_adderSum(a, b, o) {
            var x = o.s = a + b, bv = x - a, av = x - bv;
            o.t = a - av + (b - bv);
        }
        var π = Math.PI, τ = 2 * π, halfπ = π / 2, ε = 1e-6, ε2 = ε * ε, d3_radians = π / 180, d3_degrees = 180 / π;
        function d3_sgn(x) {
            return x > 0 ? 1 : x < 0 ? -1 : 0;
        }
        function d3_cross2d(a, b, c) {
            return (b[0] - a[0]) * (c[1] - a[1]) - (b[1] - a[1]) * (c[0] - a[0]);
        }
        function d3_acos(x) {
            return x > 1 ? 0 : x < -1 ? π : Math.acos(x);
        }
        function d3_asin(x) {
            return x > 1 ? halfπ : x < -1 ? -halfπ : Math.asin(x);
        }
        function d3_sinh(x) {
            return ((x = Math.exp(x)) - 1 / x) / 2;
        }
        function d3_cosh(x) {
            return ((x = Math.exp(x)) + 1 / x) / 2;
        }
        function d3_tanh(x) {
            return ((x = Math.exp(2 * x)) - 1) / (x + 1);
        }
        function d3_haversin(x) {
            return (x = Math.sin(x / 2)) * x;
        }
        d3.geo.area = function(object) {
            d3_geo_areaSum = 0;
            d3.geo.stream(object, d3_geo_area);
            return d3_geo_areaSum;
        };
        var d3_geo_areaSum, d3_geo_areaRingSum = new d3_adder();
        var d3_geo_area = {
            sphere: function() {
                d3_geo_areaSum += 4 * π;
            },
            point: d3_noop,
            lineStart: d3_noop,
            lineEnd: d3_noop,
            polygonStart: function() {
                d3_geo_areaRingSum.reset();
                d3_geo_area.lineStart = d3_geo_areaRingStart;
            },
            polygonEnd: function() {
                var area = 2 * d3_geo_areaRingSum;
                d3_geo_areaSum += area < 0 ? 4 * π + area : area;
                d3_geo_area.lineStart = d3_geo_area.lineEnd = d3_geo_area.point = d3_noop;
            }
        };
        function d3_geo_areaRingStart() {
            var λ00, φ00, λ0, cosφ0, sinφ0;
            d3_geo_area.point = function(λ, φ) {
                d3_geo_area.point = nextPoint;
                λ0 = (λ00 = λ) * d3_radians, cosφ0 = Math.cos(φ = (φ00 = φ) * d3_radians / 2 + π / 4), 
                sinφ0 = Math.sin(φ);
            };
            function nextPoint(λ, φ) {
                λ *= d3_radians;
                φ = φ * d3_radians / 2 + π / 4;
                var dλ = λ - λ0, sdλ = dλ >= 0 ? 1 : -1, adλ = sdλ * dλ, cosφ = Math.cos(φ), sinφ = Math.sin(φ), k = sinφ0 * sinφ, u = cosφ0 * cosφ + k * Math.cos(adλ), v = k * sdλ * Math.sin(adλ);
                d3_geo_areaRingSum.add(Math.atan2(v, u));
                λ0 = λ, cosφ0 = cosφ, sinφ0 = sinφ;
            }
            d3_geo_area.lineEnd = function() {
                nextPoint(λ00, φ00);
            };
        }
        function d3_geo_cartesian(spherical) {
            var λ = spherical[0], φ = spherical[1], cosφ = Math.cos(φ);
            return [ cosφ * Math.cos(λ), cosφ * Math.sin(λ), Math.sin(φ) ];
        }
        function d3_geo_cartesianDot(a, b) {
            return a[0] * b[0] + a[1] * b[1] + a[2] * b[2];
        }
        function d3_geo_cartesianCross(a, b) {
            return [ a[1] * b[2] - a[2] * b[1], a[2] * b[0] - a[0] * b[2], a[0] * b[1] - a[1] * b[0] ];
        }
        function d3_geo_cartesianAdd(a, b) {
            a[0] += b[0];
            a[1] += b[1];
            a[2] += b[2];
        }
        function d3_geo_cartesianScale(vector, k) {
            return [ vector[0] * k, vector[1] * k, vector[2] * k ];
        }
        function d3_geo_cartesianNormalize(d) {
            var l = Math.sqrt(d[0] * d[0] + d[1] * d[1] + d[2] * d[2]);
            d[0] /= l;
            d[1] /= l;
            d[2] /= l;
        }
        function d3_geo_spherical(cartesian) {
            return [ Math.atan2(cartesian[1], cartesian[0]), d3_asin(cartesian[2]) ];
        }
        function d3_geo_sphericalEqual(a, b) {
            return abs(a[0] - b[0]) < ε && abs(a[1] - b[1]) < ε;
        }
        d3.geo.bounds = function() {
            var λ0, φ0, λ1, φ1, λ_, λ__, φ__, p0, dλSum, ranges, range;
            var bound = {
                point: point,
                lineStart: lineStart,
                lineEnd: lineEnd,
                polygonStart: function() {
                    bound.point = ringPoint;
                    bound.lineStart = ringStart;
                    bound.lineEnd = ringEnd;
                    dλSum = 0;
                    d3_geo_area.polygonStart();
                },
                polygonEnd: function() {
                    d3_geo_area.polygonEnd();
                    bound.point = point;
                    bound.lineStart = lineStart;
                    bound.lineEnd = lineEnd;
                    if (d3_geo_areaRingSum < 0) λ0 = -(λ1 = 180), φ0 = -(φ1 = 90); else if (dλSum > ε) φ1 = 90; else if (dλSum < -ε) φ0 = -90;
                    range[0] = λ0, range[1] = λ1;
                }
            };
            function point(λ, φ) {
                ranges.push(range = [ λ0 = λ, λ1 = λ ]);
                if (φ < φ0) φ0 = φ;
                if (φ > φ1) φ1 = φ;
            }
            function linePoint(λ, φ) {
                var p = d3_geo_cartesian([ λ * d3_radians, φ * d3_radians ]);
                if (p0) {
                    var normal = d3_geo_cartesianCross(p0, p), equatorial = [ normal[1], -normal[0], 0 ], inflection = d3_geo_cartesianCross(equatorial, normal);
                    d3_geo_cartesianNormalize(inflection);
                    inflection = d3_geo_spherical(inflection);
                    var dλ = λ - λ_, s = dλ > 0 ? 1 : -1, λi = inflection[0] * d3_degrees * s, antimeridian = abs(dλ) > 180;
                    if (antimeridian ^ (s * λ_ < λi && λi < s * λ)) {
                        var φi = inflection[1] * d3_degrees;
                        if (φi > φ1) φ1 = φi;
                    } else if (λi = (λi + 360) % 360 - 180, antimeridian ^ (s * λ_ < λi && λi < s * λ)) {
                        var φi = -inflection[1] * d3_degrees;
                        if (φi < φ0) φ0 = φi;
                    } else {
                        if (φ < φ0) φ0 = φ;
                        if (φ > φ1) φ1 = φ;
                    }
                    if (antimeridian) {
                        if (λ < λ_) {
                            if (angle(λ0, λ) > angle(λ0, λ1)) λ1 = λ;
                        } else {
                            if (angle(λ, λ1) > angle(λ0, λ1)) λ0 = λ;
                        }
                    } else {
                        if (λ1 >= λ0) {
                            if (λ < λ0) λ0 = λ;
                            if (λ > λ1) λ1 = λ;
                        } else {
                            if (λ > λ_) {
                                if (angle(λ0, λ) > angle(λ0, λ1)) λ1 = λ;
                            } else {
                                if (angle(λ, λ1) > angle(λ0, λ1)) λ0 = λ;
                            }
                        }
                    }
                } else {
                    point(λ, φ);
                }
                p0 = p, λ_ = λ;
            }
            function lineStart() {
                bound.point = linePoint;
            }
            function lineEnd() {
                range[0] = λ0, range[1] = λ1;
                bound.point = point;
                p0 = null;
            }
            function ringPoint(λ, φ) {
                if (p0) {
                    var dλ = λ - λ_;
                    dλSum += abs(dλ) > 180 ? dλ + (dλ > 0 ? 360 : -360) : dλ;
                } else λ__ = λ, φ__ = φ;
                d3_geo_area.point(λ, φ);
                linePoint(λ, φ);
            }
            function ringStart() {
                d3_geo_area.lineStart();
            }
            function ringEnd() {
                ringPoint(λ__, φ__);
                d3_geo_area.lineEnd();
                if (abs(dλSum) > ε) λ0 = -(λ1 = 180);
                range[0] = λ0, range[1] = λ1;
                p0 = null;
            }
            function angle(λ0, λ1) {
                return (λ1 -= λ0) < 0 ? λ1 + 360 : λ1;
            }
            function compareRanges(a, b) {
                return a[0] - b[0];
            }
            function withinRange(x, range) {
                return range[0] <= range[1] ? range[0] <= x && x <= range[1] : x < range[0] || range[1] < x;
            }
            return function(feature) {
                φ1 = λ1 = -(λ0 = φ0 = Infinity);
                ranges = [];
                d3.geo.stream(feature, bound);
                var n = ranges.length;
                if (n) {
                    ranges.sort(compareRanges);
                    for (var i = 1, a = ranges[0], b, merged = [ a ]; i < n; ++i) {
                        b = ranges[i];
                        if (withinRange(b[0], a) || withinRange(b[1], a)) {
                            if (angle(a[0], b[1]) > angle(a[0], a[1])) a[1] = b[1];
                            if (angle(b[0], a[1]) > angle(a[0], a[1])) a[0] = b[0];
                        } else {
                            merged.push(a = b);
                        }
                    }
                    var best = -Infinity, dλ;
                    for (var n = merged.length - 1, i = 0, a = merged[n], b; i <= n; a = b, ++i) {
                        b = merged[i];
                        if ((dλ = angle(a[1], b[0])) > best) best = dλ, λ0 = b[0], λ1 = a[1];
                    }
                }
                ranges = range = null;
                return λ0 === Infinity || φ0 === Infinity ? [ [ NaN, NaN ], [ NaN, NaN ] ] : [ [ λ0, φ0 ], [ λ1, φ1 ] ];
            };
        }();
        if (typeof define === "function" && define.amd) {
            define(d3);
        } else if (typeof module === "object" && module.exports) {
            module.exports = d3;
        } else {
            this.d3 = d3;
        }
    }();
})();