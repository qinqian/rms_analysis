get_args = function(x) {
    require(argparse)
    parser = ArgumentParser(description='seurat normalization')

    parser$add_argument('--dir10x', dest='data')
    parser$add_argument('--dirvel', dest='vel')
    args = parser$parse_args()

    args
}

args = get_args()
