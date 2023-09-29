#########################
# table column selector #
#########################

# basic properties
show_band = True
show_group = True
show_isolated = True
show_width = True
show_gap = True
show_gap_width = True

# band topology
show_std_B_norm = True
show_C = True

# quantum geometry tensor properties
show_C_geom_tensor = False
show_av_gxx = False
show_std_gxx = False
show_av_gxy = False
show_std_gxy = False
show_T = False
show_D = False

# column groups
basic_columns = [show_band, show_group, show_isolated, show_width, show_gap, show_gap_width]
topology_columns = [show_C, show_std_B_norm]
geometry_columns = [show_C_geom_tensor, show_av_gxx, show_std_gxx, show_av_gxy, show_std_gxy, show_T, show_D]
