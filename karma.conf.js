var webpackConfig = require('./webpack.config.js');
require('babel-polyfill');

module.exports = function(config) {
  config.set({
    basePath: '',
    frameworks: ['jasmine', 'server-side'],

    files: [
      'node_modules/babel-polyfill/dist/polyfill.js',
      'node_modules/react/umd/react.development.js',
      'node_modules/react-dom/umd/react-dom.development.js',
      'node_modules/pixi.js/dist/pixi.js',
      'node_modules/react-bootstrap/dist/react-bootstrap.js',
      'node_modules/bootstrap/dist/css/bootstrap.min.css',
      'test/SVGExportTests.js',
      'node_modules/higlass/dist/hglib.css',
      'node_modules/higlass/dist/hglib.js',
    ],

    preprocessors: {
      // add webpack as preprocessor
      'src/**/*.+(js|jsx)': ['webpack', 'sourcemap'],
      'test/**/*.+(js|jsx)': ['webpack', 'sourcemap'],
    },

    webpack: webpackConfig,

    webpackServer: {
      noInfo: true, // please don't spam the console when running in karma!
    },

    plugins: [
      'karma-webpack',
      'karma-jasmine',
      'karma-sourcemap-loader',
      'karma-chrome-launcher',
      'karma-phantomjs-launcher',
      'karma-server-side',
      'karma-verbose-reporter'
    ],

    babelPreprocessor: {
      options: {
        presets: ['airbnb'],
      },
    },
    reporters: ['verbose'],
    port: 9876,
    colors: true,
    logLevel: config.LOG_DEBUG,
    autoWatch: true,
    browsers: ['Chrome'],
    singleRun: false,
    concurrency: 1, // because we use server-side
    customLaunchers: {
      Chrome_travis_ci: {
        base: 'Chrome',
        flags: ['--no-sandbox'],
      },
    },
  });

  if (process.env.TRAVIS) {
    config.browsers = ['Chrome_travis_ci'];
  }
};
