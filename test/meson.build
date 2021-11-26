test_sources = [
  'test.cc',
]

inc = include_directories('.', '../libplaneview')

test_app = executable('test_app',
                      test_sources,
		      dependencies: gtkmm_dep,
                      include_directories: inc,
                      link_with: plane_view_lib)

test('plane-view test', test_app)
