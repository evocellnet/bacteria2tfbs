#!/usr/bin/env python

# Copyright (C) <2015> EMBL-European Bioinformatics Institute

# This program is free software: you can redistribute it and/or
# modify it under the terms of the GNU General Public License as
# published by the Free Software Foundation, either version 3 of
# the License, or (at your option) any later version.

# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
# GNU General Public License for more details.

# Neither the institution name nor the name bacteria2tfbs can
# be used to endorse or promote products derived from this
# software without prior written permission. For written
# permission, please contact <marco@ebi.ac.uk>.

# Products derived from this software may not be called
# bacteria2tfbs nor may bacteria2tfbs appear in their names
# without prior written permission of the developers.
# You should have received a copy of the GNU General Public
# License along with this program. If not, see
# <http://www.gnu.org/licenses/>.

__author__ = 'Marco Galardini'
__version__ = '0.0.1'

import Bio

def get_upstreams(sequence, upstream=-1,
        downstream=0, circular=False):
    # Get all the CDS from the sequence object
    cds = filter(lambda x: x.type=='CDS' and
                 'pseudo' not in x.qualifiers,
                 sequence.features)

    pos_start = sorted([(int(g.location.start), int(g.location.end))
        for g in cds], key=lambda x: x[1])
    pos_end = sorted(pos_start, key=lambda x: x[0])

    for f in cds:
        # Get the locus_tag
        locuses = f.qualifiers.get('locus_tag', ['None'])
        locus = locuses[0]

        if f.strand > 0:
            # Get the closest start/end of a gene
            # If a start is closer, then is an overlap
            end_list = list(filter(lambda x: x[1] < int(f.location.start),
                    pos_end))
            if len(end_list) == 0:
                end = 0
            else:
                end = end_list[-1][1]
            
            start_list = list(filter(lambda x: x[0] < int(f.location.start),
                    pos_start))
            if len(start_list) == 0:
                start = 0
            else:
                start = start_list[-1][0]

            # Corner case: a gene that overlaps
            # Should be almost impossible
            # But you never know with bacteria
            overlap_list = list(filter(lambda x: x[1] > int(f.location.start)
                                            and x[0] < int(f.location.start)
                                            and (int(x[0]) != int(f.location.start)
                                            and int(x[1]) != int(f.location.end)),
                            pos_start))

            if len(overlap_list) > 0:
                up = sequence[f.location.start :
                        f.location.start + downstream].seq
                ustart = f.location.start
                uend = f.location.start + downstream
            elif end >= start:
                up = sequence[end : f.location.start + downstream].seq
                ustart = end
                uend = f.location.start + downstream
                if len(up) > upstream + downstream and upstream > 0:
                    up = sequence[f.location.start - upstream :
                            f.location.start + downstream].seq
                    ustart = f.location.start -upstream
                    uend = f.location.start + downstream
            else:
                up = sequence[f.location.start :
                        f.location.start + downstream].seq
                ustart = f.location.start
                uend = f.location.start + downstream
        
        else:
            # Get the closest start/end of a gene
            # If a end is closer, then is an overlap
            end_list = list(filter(lambda x: x[1] > int(f.location.end),
                    pos_end))
            if len(end_list) == 0:
                end = 0
            else:
                end = end_list[0][1]
            
            start_list = list(filter(lambda x: x[0] > int(f.location.end),
                    pos_start))
            if len(start_list) == 0:
                start = 0
            else:
                start = start_list[0][0]
            
            # Corner case: a previous gene that overlaps
            # Should be almost impossible
            # But you never know with bacteria
            overlap_list = list(filter(lambda x: x[1] > int(f.location.end)
                                            and x[0] < int(f.location.end)
                                            and (int(x[0]) != int(f.location.start)
                                            and int(x[1]) != int(f.location.end)),
                            pos_start))
            
            if len(overlap_list) > 0:
                up = sequence[f.location.end - downstream :
                        f.location.end].reverse_complement().seq
                ustart = f.location.end - downstream
                uend = f.location.end
            elif start <= end:
                up = sequence[f.location.end - downstream :
                        start].reverse_complement().seq
                ustart = f.location.end - downstream
                uend = start
                if len(up) > upstream + downstream and upstream > 0:
                    up = sequence[f.location.end - downstream :
                            f.location.end + upstream].reverse_complement().seq
                    ustart = f.location.end - downstream
                    uend = f.location.end + upstream
            else:
                up = sequence[f.location.end - downstream :
                        f.location.end].reverse_complement().seq
                ustart = f.location.end - downstream
                uend = f.location.end

        yield locus, up, int(ustart), int(uend), f.strand

