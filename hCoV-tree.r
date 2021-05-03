tree <- read.tree("SARS-CoV-2/genome_sequences.nwk")
p <- ggtree(tree)
table(p$data$isTip)

p2 <- function(number, size, type) {
    p1 <- ggtree(tree, layout = type) + xlim(0, 0.02) + geom_tiplab(aes(angle = angle), size = size)
    p2 <- p1 %>% collapse(node= number) +
    geom_point2(aes(subset=(node== number)), shape = 21, size = size, fill='green')
    print(p2)
    ggsave(file = "/tmp/p2.pdf")
}
p2(208, 3, "circular")

p3 <- function(number, size, type) {
    p1 <- ggtree(tree, layout = type) + xlim(0, 0.02)
    p2 <- p1 %>% collapse(node= number) +
    geom_point2(aes(subset=(node== number)), shape = 21, size = size, fill='green')
    print(p2)
    ggsave(file = "/tmp/p3.pdf")
}
p3(208, 3, "circular")

p1 <- ggtree(tree, layout = "circular") + xlim(0, 0.1) + geom_tiplab(aes(angle = angle), size = 2)
p2 <- p1 %>% collapse(node= number) +
geom_point2(aes(subset=(node== number)), shape = 21, size = size, fill='green')
print(p2)
