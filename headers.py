clouds_path = '/home/fmlopez/'

clouds_info = [
    {
        'src': f'{clouds_path}/INAER_2011_Alcoy.xyz',
        'dst': f'{clouds_path}/small_AlcoyH.xyz',
        'Npoints': 2772832,
        'minx': 715244.96,
        'maxx': 716057.75,
        'miny': 4286623.63,
        'maxy': 4287447.70
    },
    {
        'src': f'{clouds_path}/INAER_2011_Alcoy_Core.xyz',
        'dst': f'{clouds_path}/AlcoyH.xyz',
        'Npoints': 20380212,
        'minx': 714947.98,
        'maxx': 716361.06,
        'miny': 4286501.93,
        'maxy': 4288406.23
    },
    {
        'src': f'{clouds_path}/BABCOCK_2017_Arzua_3B.xyz',
        'dst': f'{clouds_path}/ArzuaH.xyz',
        'Npoints': 40706503,
        'minx': 568000.00,
        'maxx': 568999.99,
        'miny': 4752320.00,
        'maxy': 4753319.99
    },
    {
        'src': f'{clouds_path}/V19_group1_densified_point_cloud.xyz',
        'dst': f'{clouds_path}/BrionUH.xyz',
        'Npoints': 48024480,
        'minx': 526955.908,
        'maxx': 527686.445,
        'miny': 4742586.025,
        'maxy': 4743124.373
    },
    {
        'src': f'{clouds_path}/V21_group1_densified_point_cloud.xyz',
        'dst': f'{clouds_path}/BrionFH.xyz',
        'Npoints': 42384876,
        'minx': 526964.093,
        'maxx': 527664.647,
        'miny': 4742610.292,
        'maxy': 4743115.738
    },
]

for cloud in clouds_info:

    src = cloud['src']
    dst = cloud['dst']
    Npoints = cloud['Npoints']
    minx = cloud['minx']
    maxx = cloud['maxx']
    miny = cloud['miny']
    maxy = cloud['maxy']

    # insert a line at the beginning of the file
    with open(dst, 'w') as file_dst:
        file_dst.write(f'{Npoints} {minx} {maxx} {miny} {maxy}\n')
        # copy the rest of the file
        with open(src, 'r') as file_src:
            for line in file_src:
                file_dst.write(line)

    print(f'Header added to file {dst}')