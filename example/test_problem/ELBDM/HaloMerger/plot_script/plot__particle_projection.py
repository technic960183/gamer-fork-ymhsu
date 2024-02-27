import argparse
import sys
import yt

# load the command-line parameters
parser = argparse.ArgumentParser( description='Projection of the particles' )

parser.add_argument( '-p', action='store', required=False, type=str, dest='prefix',
                     help='path prefix [%(default)s]', default='../' )
parser.add_argument( '-s', action='store', required=True,  type=int, dest='idx_start',
                     help='first data index' )
parser.add_argument( '-e', action='store', required=True,  type=int, dest='idx_end',
                     help='last data index' )
parser.add_argument( '-d', action='store', required=False, type=int, dest='didx',
                     help='delta data index [%(default)d]', default=1 )

args=parser.parse_args()

# take note
print( '\nCommand-line arguments:' )
print( '-------------------------------------------------------------------' )
print( ' '.join(map(str, sys.argv)) )
print( '-------------------------------------------------------------------\n' )


idx_start   = args.idx_start
idx_end     = args.idx_end
didx        = args.didx
prefix      = args.prefix

colormap    = 'algae'
dpi         = 150

yt.enable_parallelism()

ts = yt.DatasetSeries( [ prefix+'/Data_%06d'%idx for idx in range(idx_start, idx_end+1, didx) ] )
for ds in ts.piter():
   for center_mode in ['c', 'm']:
      for direction in ['y', 'z']:
         for field in ['particle_mass']:

            # decide the center
            if center_mode == 'm':
               center = ds.all_data().quantities.max_location('ParDens')[1:]
            elif center_mode == 'c':
               center = 'c'

            # decide the direction
            if direction == 'z':
               p_par = yt.ParticlePlot( ds, 'particle_position_x', 'particle_position_y', field, center=center )
            elif direction == 'y':
               p_par = yt.ParticlePlot( ds, 'particle_position_z', 'particle_position_x', field, center=center )

            # setting for the figure
            p_par.set_axes_unit( 'kpc' )
            p_par.set_unit( field, 'Msun'        )
            p_par.set_zlim( field, 1.0e+2, 1.0e7 )
            p_par.set_cmap( field, colormap      )
            p_par.set_background_color( field )
            p_par.annotate_timestamp( time_unit='Gyr', corner='upper_right' )

            # zoom in
            if center_mode == 'm':
               p_par.zoom(4)

            # save the figure
            p_par.save( '%s_%s'%(ds, center_mode), mpl_kwargs={'dpi':dpi} )
