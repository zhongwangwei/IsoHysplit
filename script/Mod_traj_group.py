from __future__ import division, print_function

import os
import numpy as np

from Mod_traj_util import *

class HyGroup(object):
    """
    Class for initializing ``HyPath`` container object.

    :superclass: for ``TrajectoryGroup`` and ``ClusterGroup``.

    """

    def __init__(self, trajectories):
        """
        Initialize ``HyGroup``.

        Parameters
        ----------
        trajectories : list of ``Trajectory`` instances
            ``Trajectory`` instances that belong to the group

        """
        self.trajectories = trajectories
        self.trajcount = len(trajectories)
        self.trajids = [traj.trajid for traj in self.trajectories]

    def __add__(self, other):
        """
        Add ``HyGroup`` instances.

        Create ``HyGroup`` from the union of two sets of ``Trajectory``
        instances.

        Parameters
        ----------
        other : ``HyGroup`` subclass instance

        Returns
        -------
        newgroup : list
            A list of the unique trajectories from the two groups.
            Used to make new ``TrajectoryGroup`` instance.

        """
        set0 = set(self.trajectories)
        set1 = set(other.trajectories)

        newgroup = list(set0 | set1)

        return newgroup

    def __sub__(self, other):
        """
        Subtract ``HyGroup`` instances.

        Create new ``HyGroup`` from the set difference of two
        sets of ``Trajectory`` instances.

        Parameters
        ----------
        other : ``HyGroup`` subclass instance

        Returns
        -------
        newgroup : list
            A list of the set difference of the trajectories.
            Has has all the elements of ``self`` with the
            trajectories of ``other`` removed.  Used to
            make new ``TrajectoryGroup`` instance.

        """
        set0 = set(self.trajectories)
        set1 = set(other.trajectories)

        newgroup = list(set0 - set1)

        return newgroup

    def make_infile(self, infile_dir, use_clippedpath=True):
        """
        Write member ``HyPath`` file paths to INFILE.

        INFILE is used by ``HYSPLIT`` to perform cluster analysis.
        If a specific subset of ``HyPath`` instances is needed,
        create a new ``HyGroup`` containing only qualifying
        ``HyPath`` instances.

        Parameters
        ----------
        infile_dir : string
            The directory in which to create INFILE

        use_clippedpath : Boolean
            Default True. Write out path of clipped trajectory
            rather than original trajectory.

        """
        with open(os.path.join(infile_dir, 'INFILE'), 'w') as infile:

            for traj in self:
                skip_output = False
                if traj.multitraj:
                    try:
                        output = traj.cfullpath
                    except AttributeError:
                        print(traj.trajid, " missing clusterable file")
                        skip_output = True
                else:
                    if use_clippedpath:
                        try:
                            output = traj.cfullpath
                        except AttributeError:
                            output = traj.fullpath
                    else:
                        output = traj.fullpath

                if not skip_output:
                    output = output.replace('\\', '/')
                    infile.writelines(output + '\n')
                    infile.flush()

class TrajectoryGroup(HyGroup):
    """
    Class for processing and plotting multiple ``Trajectory`` instances.

    :subclass: of ``HyGroup``.

    """

    def __init__(self, trajectories):
        """
        Initialize ``TrajectoryGroup`` object.

        Parameters
        ----------
        trajectories : list of ``Trajectory`` instances
            ``Trajectory`` instances that belong in the group.

        """
        HyGroup.__init__(self, trajectories)

    def __getitem__(self, index):
        """
        Get ``Trajectory`` or ``TrajectoryGroup``.

        Parameters
        ----------
        index : int or slice

        Returns
        -------
        ``Trajectory`` or ``TrajectoryGroup`` depending if indexed
        or sliced.  Won't return a ``Cluster`` because those are
        specially defined.

        """
        newthing = self.trajectories[index]

        if isinstance(newthing, list):
            newthing = TrajectoryGroup(newthing)

        return newthing

    def __add__(self, other):
        """
        Add a ``HyGroup`` to this ``TrajectoryGroup`` instance.

        Parameters
        ----------
        other : ``HyGroup``
            Another ``TrajectoryGroup`` or ``Cluster``.  May or may not
            contain some of the same ``Trajectory`` instances

        Returns
        -------
        A new ``TrajectoryGroup`` containing the union of the sets
        of ``Trajectory`` instances.

        """
        return TrajectoryGroup(HyGroup.__add__(self, other))

    def __sub__(self, other):
        """
        Subtract a ``HyGroup`` from this ``TrajectoryGroup`` instance.

        Parameters
        ----------
        other : ``HyGroup``
            Another ``TrajectoryGroup`` or ``Cluster``

        Returns
        -------
        A new ``TrajectoryGroup`` containing the set difference between
        the sets of ``Trajectory`` instances.

        """
        return TrajectoryGroup(HyGroup.__sub__(self, other))

    def pop(self, ind=-1, trajid=None):
        """
        Remove Trajectory object(s) from self.

        Shortcut to self.trajectories.pop() that updates the
        self.trajcount and the list of trajids.

        Parameters
        ----------
        ind : int
            The positional argument of the ``Trajectory``
            to remove.
        trajid : string or list of strings
            The identifier(s) of the ``Trajectory`` object(s)
            to remove from ``self``.  Overrides ``ind`` if not None.

        Returns
        -------
        popped : ``Trajectory`` or ``TrajectoryGroup``
            A``Trajectory`` or ``TrajectoryGroup`` consisting of the
            trajectory or trajectories indicated by ``ind`` or ``trajid``.

        """
        if trajid is not None:
            try:
                to_pop = [self.trajids.index(trajid)]
            except ValueError:
                to_pop = [self.trajids.index(t) for t in trajid
                          if t in self.trajids]
                if len(to_pop) == 0:
                    raise ValueError('TrajIDs not in list of self.trajids')

            to_pop.sort()
            popped = []
            for p in to_pop[::-1]:
                popped.append(self.trajectories.pop(p))
                self.trajids.pop(p)
            self.trajcount = len(self.trajectories)

            if len(popped) == 1:
                popped = popped[0]
            else:
                popped = TrajectoryGroup(popped)
        else:
            popped = self.trajectories.pop(ind)
            self.trajids.pop(ind)
            self.trajcount = len(self.trajectories)

        return popped

    def append(self, traj):
        """
        Add a ``Trajectory`` to the ``self``.

        Parameters
        ----------
        traj : ``Trajectory`` instance
            The ``Trajectory`` to add to the end of ``self``.

        """
        if hasattr(traj, 'trajid'):
            self.trajectories.append(traj)
            self.trajids.append(traj.trajid)
            self.trajcount = len(self.trajectories)

def make_trajectorygroup(signature):
    """
    Initialize ``Trajectory`` instances from HYSPLIT trajectory data files.

    Parameters
    ----------
    signature : string or iterable of strings
        Signature shared by a group of HYSPLIT simulation files from one or
        multiple model runs (if multiple, must contain same output variables).
        This is a Bash-style signature, not a real expression.  The `*` char is
        a wildcard. Can include an absolute or relative path, or no path to
        target the current directory. If ``signature`` is not a string, it is
        assumed to be an iterable containing path(s) to specifc HYSPLIT files.

    Returns
    -------
    trajectories : ``TrajectoryGroup``
        ``TrajectoryGroup`` object containing ``Trajectory`` instances created
        from all simulation files matching ``signature``.

    """
    # Get list of hysplit files matching signature
    if 'basestring' not in globals():  # Python 2/3 compat
        basestring = str
    if not isinstance(signature, basestring):
        hyfiles = [os.path.split(hyfile)[-1] for hyfile in signature]
        folder, _ = os.path.split(signature[0])
    else:
        hyfiles = hysplit_filelister(signature)
        folder, _ = os.path.split(signature)

    orig_dir = os.getcwd()
    trajectories = []

    # Wrap in try .. finally to ensure current working directory restored
    try:
        os.chdir(folder)

        # Sort list of hysplit files by the datestring at the end
        # Will also sort in ascending altitude within each same datestring
        hyfiles.sort(key=lambda x: x[-8:])

        clipdir = os.path.join(folder, 'clippedtraj')
        if not os.path.isdir(clipdir):
            clipdir = None

        # Load in the hysplit file data
        # Get lists of the datestrings and filenames of the hysplit files
        for hyfile in hyfiles:

            data, path, head, datetime, multitraj = load_hysplitfile(hyfile)

            if multitraj:
                # Initialize trajectory objects
                for d, p, dt in zip(data, path, datetime):

                    # Get rid of parcel number in d
                    # Get rid of parcel #, lat, lon, altitude in head
                    trajectories.append(Trajectory(d, p, dt, head,
                                                   folder, hyfile, clipdir,
                                                   multitraj))

            else:
                trajectories.append(Trajectory(data, path, datetime, head,
                                               folder, hyfile, clipdir,
                                               multitraj))

        # initialize trajectory group
        trajectories = TrajectoryGroup(trajectories)
    finally:
        # Restore current working directory no matter what
        os.chdir(orig_dir)

    return trajectories

