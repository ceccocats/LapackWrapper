%w(colorize rake fileutils).each do |gem|
  begin
    require gem
  rescue LoadError
    warn "Install the #{gem} gem:\n $ (sudo) gem install #{gem}".magenta
    exit 1
  end
end

task :default => [:build]

task :mkl, [:year, :bits] do |t, args|
  args.with_defaults(:year => "2017", :bits => "x64" )
  sh "'C:/Program Files (x86)/IntelSWTools/compilers_and_libraries/windows/bin/compilervars.bat' -arch #{args.bits} vs#{args.year}shell"
end

TESTS = [
  "test1-small-factorization",
  "test2-Timing",
  "test3-BandedMatrix",
  "test4-BFGS",
  "test5-BLOCKTRID",
  "test6-EIGS"
]

desc "run tests"
task :run do
  TESTS.each do |cmd|
    sh "./bin/#{cmd}"
  end
end

task :run_win do
  TESTS.each do |cmd|
    sh "bin\\Release\\#{cmd}.exe"
  end
end

desc "build lib"
task :build  do
  sh "make"
end

def ChangeOnFile( file, text_to_replace, text_to_put_in_place )
  text = File.read file
  File.open(file, 'w+'){|f| f << text.gsub(text_to_replace, text_to_put_in_place)}
end

desc "compile for Visual Studio [default year=2017 bits=x64]"
task :build_win, [:year, :bits, :lapack] do |t, args|
  args.with_defaults(
    :year   => "2017",
    :bits   => "x64",
    :lapack => "LAPACK_WRAPPER_USE_OPENBLAS"
    #:lapack => "LAPACK_WRAPPER_USE_LAPACK",
    #:lapack => "LAPACK_WRAPPER_USE_MKL"
  )

  cmd = "set path=%path%;lib3rd\\lib;lib3rd\\dll;"

  FileUtils.rm_f 'src/lapack_wrapper/lapack_wrapper_config.hh'
  FileUtils.cp   'src/lapack_wrapper/lapack_wrapper_config.hh.tmpl', 'src/lapack_wrapper/lapack_wrapper_config.hh'

  ChangeOnFile(
    'src/lapack_wrapper/lapack_wrapper_config.hh',
    '@@LAPACK_WRAPPER_USE@@',
    "#define #{args.lapack} 1"
  )
  ChangeOnFile(
    'src/lapack_wrapper/lapack_wrapper_config.hh',
    '@@LAPACK_WRAPPER_THREAD@@',
    "#define #{args.thread} 1"
  )
  ChangeOnFile(
    'src/lapack_wrapper/lapack_wrapper_config.hh',
    '@@LAPACK_WRAPPER_NOSYSTEM_OPENBLAS@@',
    "#define LAPACK_WRAPPER_DO_NOT_USE_SYSTEM_OPENBLAS 1"
  )

  dir = "vs_#{args.year}_#{args.bits}"

  FileUtils.rm_rf   dir
  FileUtils.mkdir_p dir
  FileUtils.cd      dir

  # do not build executable
  #tmp = " -DBITS=#{args.bits} -DYEAR=#{args.year} " + ' -DBUILD_EXECUTABLE=1 -DCMAKE_INSTALL_PREFIX:PATH=..\lib ..'
  tmp = " -DBITS=#{args.bits} -DYEAR=#{args.year} " +
        ' -DCMAKE_INSTALL_PREFIX:PATH=..\lib ' +
        ' -DBUILD_EXECUTABLE:VAR=true ..'

  win32_64 = ''
  case args.bits
  when /x64/
    win32_64 = ' Win64'
  end

  case args.year
  when "2010"
    sh 'cmake -G "Visual Studio 10 2010' + win32_64 +'" ' + tmp
  when "2012"
    sh 'cmake -G "Visual Studio 11 2012' + win32_64 +'" ' + tmp
  when "2013"
    sh 'cmake -G "Visual Studio 12 2013' + win32_64 +'" ' + tmp
  when "2015"
    sh 'cmake -G "Visual Studio 14 2015' + win32_64 +'" ' + tmp
  when "2017"
    sh 'cmake -G "Visual Studio 15 2017' + win32_64 +'" ' + tmp
  else
    puts "Visual Studio year #{year} not supported!\n";
  end

  sh 'cmake --build . --config Release  --target install'
  FileUtils.mkdir_p "../lib/lib"
  FileUtils.mkdir_p "../lib/bin"
  FileUtils.mkdir_p "../lib/dll"
  FileUtils.mkdir_p "../lib/include"
  FileUtils.cp "lib/bin/HSL_#{args.bits}.dll",                   "../lib/bin/libHSL_#{args.bits}.dll"
  FileUtils.cp "lib/lib/HSL_#{args.bits}.lib",                   "../lib/dll/libHSL_#{args.bits}.lib"
  FileUtils.cp "lib/bin/lapack_wrapper_#{args.bits}.dll",        "../lib/bin/liblapack_wrapper_#{args.bits}.dll"
  FileUtils.cp "lib/lib/lapack_wrapper_#{args.bits}.lib",        "../lib/dll/liblapack_wrapper_#{args.bits}.lib"
  FileUtils.cp "lib/lib/lapack_wrapper_#{args.bits}_static.lib", "../lib/lib/liblapack_wrapper_#{args.bits}_static.lib"

  FileUtils.cp_r "lib/include", "../lib/"

  sh 'cmake --build . --config Debug --target install'
  FileUtils.cp "lib/bin/HSL_#{args.bits}.dll",            "../lib/bin/libHSL_#{args.bits}_debug.dll"
  FileUtils.cp "lib/lib/HSL_#{args.bits}.lib",            "../lib/dll/libHSL_#{args.bits}.lib"
  FileUtils.cp "lib/bin/lapack_wrapper_#{args.bits}.dll", "../lib/bin/liblapack_wrapper_#{args.bits}_debug.dll"
  FileUtils.cp "lib/lib/lapack_wrapper_#{args.bits}.lib", "../lib/dll/liblapack_wrapper_#{args.bits}_debug.lib"
  FileUtils.cp "lib/lib/lapack_wrapper_#{args.bits}_static.lib", "../lib/lib/liblapack_wrapper_#{args.bits}_static_debug.lib"

  FileUtils.cd '..'

end

desc 'compile for OSX [default lapack="LAPACK_WRAPPER_USE_ACCELERATE"]'
task :build_osx, [:lapack] do |t, args|
  args.with_defaults(
    :lapack => "LAPACK_WRAPPER_USE_ACCELERATE"
    #:lapack => "LAPACK_WRAPPER_USE_LAPACK",
    #:lapack => "LAPACK_WRAPPER_USE_MKL"
  )

  FileUtils.rm_f 'src/lapack_wrapper/lapack_wrapper_config.hh'
  FileUtils.cp   'src/lapack_wrapper/lapack_wrapper_config.hh.tmpl', 'src/lapack_wrapper/lapack_wrapper_config.hh'

  ChangeOnFile(
    'src/lapack_wrapper/lapack_wrapper_config.hh',
    '@@LAPACK_WRAPPER_USE@@',
    "#define #{args.lapack} 1"
  )
  ChangeOnFile(
    'src/lapack_wrapper/lapack_wrapper_config.hh',
    '@@LAPACK_WRAPPER_THREAD@@',
    "#define #{args.thread} 1"
  )
  ChangeOnFile(
    'src/lapack_wrapper/lapack_wrapper_config.hh',
    '@@LAPACK_WRAPPER_NOSYSTEM_OPENBLAS@@',
    "// #define LAPACK_WRAPPER_DO_NOT_USE_SYSTEM_OPENBLAS 1"
  )

  dir = "build"

  FileUtils.rm_rf   dir
  FileUtils.mkdir_p dir
  FileUtils.cd      dir

  # do not build executable
  sh 'cmake ..'
  sh 'cmake --build . --config Release --target install'
  FileUtils.mkdir_p "../lib"
  FileUtils.cp_r    './lib/lib',     '../lib/'
  #FileUtils.cp_r    './lib/bin',     '../lib/'
  FileUtils.cp_r    './lib/include', '../lib/'
  FileUtils.cd '..'

end

desc 'compile for LINUX [default lapack="LAPACK_WRAPPER_USE_OPENBLAS"]'
task :build_linux, [:lapack] do |t, args|
  args.with_defaults(
    :lapack => "LAPACK_WRAPPER_USE_OPENBLAS"
    #:lapack => "LAPACK_WRAPPER_USE_LAPACK",
    #:lapack => "LAPACK_WRAPPER_USE_MKL"
  )

  FileUtils.rm_f 'src/lapack_wrapper/lapack_wrapper_config.hh'
  FileUtils.cp   'src/lapack_wrapper/lapack_wrapper_config.hh.tmpl', 'src/lapack_wrapper/lapack_wrapper_config.hh'

  ChangeOnFile(
    'src/lapack_wrapper/lapack_wrapper_config.hh',
    '@@LAPACK_WRAPPER_USE@@',
    "#define #{args.lapack} 1"
  )
  ChangeOnFile(
    'src/lapack_wrapper/lapack_wrapper_config.hh',
    '@@LAPACK_WRAPPER_THREAD@@',
    "#define #{args.thread} 1"
  )
  ChangeOnFile(
    'src/lapack_wrapper/lapack_wrapper_config.hh',
    '@@LAPACK_WRAPPER_NOSYSTEM_OPENBLAS@@',
    "// #define LAPACK_WRAPPER_DO_NOT_USE_SYSTEM_OPENBLAS 1"
  )

  dir = "build"

  FileUtils.rm_rf   dir
  FileUtils.mkdir_p dir
  FileUtils.cd      dir

  # do not build executable
  sh 'cmake ..'
  sh 'cmake --build . --config Release --target install'
  FileUtils.mkdir_p "../lib"
  FileUtils.mkdir_p "../lib/lib"
  FileUtils.mkdir_p "../lib/bin"
  FileUtils.cp_r    './lib/lib',     '../lib/'
  FileUtils.cp_r    './lib/bin',     '../lib/'
  FileUtils.cp_r    './lib/include', '../lib/'
  FileUtils.cd '..'

end
