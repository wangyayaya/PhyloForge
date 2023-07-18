#!/usr/bin/env python
# -*- coding:utf-8 -*-

import argparse
import os.path
import sys
import configparser

'''sys.path.append(os.getcwd())

import est
import est.mode as mode
'''
import mode


def get_parser():
    parser = argparse.ArgumentParser(description=f"To get started: {__file__} -c run.config")
    parser.add_argument("-c", "--config", type=str, help=f"run main program, You can obtain configfile by running: "
                                                         f"{__file__} -g >run.cfg, or refer README.md")
    parser.add_argument("-grc", "--get_run_cfg", action="store_true", help="get the main config file")
    parser.add_argument("-gsc", "--get_soft_cfg", action="store_true", help="get the software configuration file")
    parser.add_argument("-sc", "--soft_cfg", help="perform software configuration")
    parser.add_argument("-v", "--version", action="version", version="est 测试版")
    args = parser.parse_args()

    if args.get_run_cfg:
        cfg = os.path.join(est.__path__[0], "cfg/run.config")
        with open(cfg) as f:
            for line in f.readlines():
                print(line, end='')
        sys.exit()

    elif args.get_soft_cfg:
        cfg = os.path.join(est.__path__[0], "cfg/software.ini")
        with open(cfg) as f:
            for line in f.readlines():
                print(line, end='')
        sys.exit()

    elif args.config:
        opt_cfg = args.config

    elif args.soft_cfg:
        conf = configparser.ConfigParser()
        conf.read(args.soft_cfg)
        ob_path = os.path.split(__file__)[0]
        cfg = os.path.join(ob_path, "cfg", "software.ini")
        if conf.sections()[0] == 'software':
            conf.write(open(os.path.join(cfg), 'w'))
            print('Software config file has been modified')
        else:
            print('Software config file no change')
        sys.exit()

    else:
        print(parser.format_usage())
        sys.exit()

    return opt_cfg


def main():
    mode.MODE().run()


if __name__ == '__main__':
    main()
