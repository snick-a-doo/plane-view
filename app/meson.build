plane_view_sources = [
  'main.cc',
]

inc = include_directories('../libplaneview')

plane_view_app = executable('plane-view',
                            plane_view_sources,
                            include_directories : inc,
                            dependencies : gtkmm_dep,
                            link_with : [plane_view_lib],
                            install : true)
