NR==1{
        next
} {
        printf("%s",$0)
} END {
        printf("\n")
}
