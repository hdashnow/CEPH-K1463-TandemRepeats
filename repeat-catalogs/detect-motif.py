

def find_motif(sequence, motif):
    """Find the edges of a given motif in the sequence. Report motif count.
    Arguments:
        sequence (str): DNA sequence to be searched
        motif (str): DNA motif be queried
    Returns:
        Coordinates of motif in sequence
        Exact motif count
    """
    window_size = 3
    window_counts = []
    coords = []
    prev_count = 0
    for i in range(0, len(sequence) - window_size, window_size):
        #print(i)
        window_seq = sequence[i:i + window_size]
        #print(window_seq)
        window_count = window_seq.count(motif)
        if i > 0:
            prev_count = window_counts[-1]
        if window_count != prev_count:
            coords.append(i)
        window_counts.append(window_count)
    print(window_counts)

    
    motif_count = sum(window_counts)
    return coords, motif_count

# Check if this is in fact the predominant motif. If others, find them, work out edges, and count them.

print(find_motif('AAACAGCAGCAGCAGCAGCAGCAGAAAAAAAAAAAACAGCAGCAGCAGCAGTTTTTTTTTTTTTTT', 'CAG'))