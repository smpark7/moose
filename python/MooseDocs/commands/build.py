#pylint: disable=missing-docstring
####################################################################################################
#                                    DO NOT MODIFY THIS HEADER                                     #
#                   MOOSE - Multiphysics Object Oriented Simulation Environment                    #
#                                                                                                  #
#                              (c) 2010 Battelle Energy Alliance, LLC                              #
#                                       ALL RIGHTS RESERVED                                        #
#                                                                                                  #
#                            Prepared by Battelle Energy Alliance, LLC                             #
#                               Under Contract No. DE-AC07-05ID14517                               #
#                               With the U. S. Department of Energy                                #
#                                                                                                  #
#                               See COPYRIGHT for full restrictions                                #
####################################################################################################
#pylint: enable=missing-docstring

import os
import math
import multiprocessing
import shutil
from distutils.dir_util import copy_tree
import logging

import livereload

import MooseDocs
from MooseDocsNode import MooseDocsNode
from MarkdownNode import MarkdownNode

LOG = logging.getLogger(__name__)

def build_options(parser):
    """
    Command-line options for build command.
    """
    parser.add_argument('--config-file', type=str, default='website.yml',
                        help="The configuration file to use for building the documentation using "
                             "MOOSE. (Default: %(default)s)")
    parser.add_argument('--num-threads', '-j', type=int, default=multiprocessing.cpu_count(),
                        help="Specify the number of threads to build pages with.")
    parser.add_argument('--template', type=str, default='website.html',
                        help="The template html file to utilize (Default: %(default)s).")

    parser.add_argument('--host', default='127.0.0.1', type=str,
                        help="The local host location for live web server (default: %(default)s).")
    parser.add_argument('--port', default='8000', type=str,
                        help="The local host port for live web server (default: %(default)s).")
    parser.add_argument('--site-dir', type=str, default=os.path.join(MooseDocs.ROOT_DIR, 'site'),
                        help="The location to build the website content (Default: %(default)s).")
    parser.add_argument('--serve', action='store_true',
                        help="Serve the presentation with live reloading, the 'site_dir' is "
                             "ignored for this case.")
    parser.add_argument('--no-livereload', action='store_true',
                        help="When --serve is used this flag disables the live reloading.")


def make_tree(directory, node, site_dir, parser):
    """
    Create the tree structure of NavigationNode/MarkdownNode objects
    """
    for p in sorted(os.listdir(directory)):
        child = None
        path = os.path.join(directory, p)
        if p in ['index.md', 'index.html']:
            continue

        if os.path.isfile(path) and (path.endswith('.md')):
            name = os.path.basename(path)[:-3]
            MarkdownNode(name=name, parent=node, markdown=path, site_dir=site_dir, parser=parser)

        elif os.path.isdir(path) and (p not in ['.', '..']):
            name = os.path.basename(path)
            md = os.path.join(path, 'index.md')
            #pylint: disable=redefined-variable-type
            if os.path.exists(md):
                child = MarkdownNode(name=name, parent=node, markdown=md, site_dir=site_dir,
                                     parser=parser)
            else:
                child = MooseDocsNode(name=name, parent=node, site_dir=site_dir)
            #pylint: enable=redefined-variable-type
            make_tree(path, child, site_dir, parser)

def flat(node):
    """
    Create a flat list of pages for parsing and generation.

    Args:
      node[NavigationNode]: The root node to flatten from
    """
    for child in node:
        if isinstance(child, MarkdownNode):
            yield child
        for c in flat(child):
            yield c

class Builder(object):
    """
    Object for building
    """
    def __init__(self, parser, site_dir):

        self._site_dir = site_dir
        md_file = os.path.join(os.getcwd(), 'content', 'index.md')
        self._root = MarkdownNode(name=str(), markdown=md_file, parser=parser,
                                  site_dir=self._site_dir)
        make_tree(os.path.dirname(md_file), self._root, self._site_dir, parser)
        self._pages = [self._root] + list(flat(self._root))

    def __iter__(self):
        """
        Allow direct iteration over pages contained in this object.
        """
        return self._pages.__iter__()

    def build(self, num_threads=multiprocessing.cpu_count()):
        """
        Build all the pages in parallel.
        """

        def make_chunks(local, num):
            """Divides objects into equal size containers for parallel execution."""
            num = int(math.ceil(len(local)/float(num)))
            for i in range(0, len(local), num):
                yield local[i:i + num]

        def build_pages(pages, lock):
            """Loops through given pages and call build method."""
            for page in pages:
                page.build(lock)

        jobs = []
        lock = multiprocessing.Lock()
        for chunk in make_chunks(self._pages, num_threads):
            p = multiprocessing.Process(target=build_pages, args=(chunk, lock))
            p.start()
            jobs.append(p)

        for job in jobs:
            job.join()

        self.copyFiles()

    def copyFiles(self):
        """
        Copy the css/js/fonts/media files for this project.
        """

        def helper(src, dst):
            """Copy helper."""
            if not os.path.exists(dst):
                os.makedirs(dst)
            if os.path.exists(src):
                copy_tree(src, dst)

        # Copy js/css/media from MOOSE and current projects
        for from_dir in [os.path.join(MooseDocs.MOOSE_DIR, 'docs'), os.getcwd()]:
            helper(os.path.join(from_dir, 'js'), os.path.join(self._site_dir, 'js'))
            helper(os.path.join(from_dir, 'css'), os.path.join(self._site_dir, 'css'))
            helper(os.path.join(from_dir, 'media'), os.path.join(self._site_dir, 'media'))

def build(config_file=None, site_dir=None, num_threads=None, no_livereload=False,
          clean=False, serve=False, host=None, port=None, **kwargs):
    """
    The main build command.
    """

    if serve:
        clean = True
        site_dir = os.path.abspath(os.path.join(MooseDocs.TEMP_DIR, 'site'))

    # Clean/create site directory
    if clean and os.path.exists(site_dir):
        LOG.info('Cleaning build directory: %s', site_dir)
        shutil.rmtree(site_dir)

    # Create the "temp" directory
    if not os.path.exists(site_dir):
        os.makedirs(site_dir)

    # Load the YAML configuration file
    config = MooseDocs.load_config(config_file, **kwargs)

    # Create the markdown parser
    parser = MooseDocs.MooseMarkdown(extensions=config.keys(), extension_configs=config)

    # Create object for storing pages to be generated
    def build_complete():
        """Builds complete documentation."""
        builder = Builder(parser, site_dir)
        builder.build(num_threads=num_threads)
        return builder
    builder = build_complete()

    # Serve
    if serve:
        # Create the live server
        server = livereload.Server()

        # Watch markdown files
        if not no_livereload:
            for page in builder:
                server.watch(page.source(), page.build)

            # Watch support directories
            server.watch(os.path.join(os.getcwd(), 'media'), builder.copyFiles)
            server.watch(os.path.join(os.getcwd(), 'css'), builder.copyFiles)
            server.watch(os.path.join(os.getcwd(), 'js'), builder.copyFiles)
            server.watch(os.path.join(os.getcwd(), 'fonts'), builder.copyFiles)

            # Watch the files and directories that require complete rebuild
            server.watch(config_file, build_complete)
            server.watch('templates', builder.build)

        # Start the server
        server.serve(root=site_dir, host=host, port=port, restart_delay=0)

    return None
