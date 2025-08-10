test_sources = [
  '../app/axis.cc',
  '../app/color.cc',
  '../app/plotter.cc',
  'test.cc',
]

inc = include_directories('.', '../app')

test_app = executable('test_app',
                      test_sources,
                      dependencies: gtkmm_dep,
                      include_directories: inc)

test('plane-view test', test_app)
