
# GET CONTACTS: Rows are participants, cols or contacts
# Conversation contacts
contactMatrixCon <- matrix(read_table(here::here("inst", "extdata", "cntCA.txt"), col_names = FALSE )$`X1`, 25, 25)  # contact matrix for physical 
# Physical contacts
contactMatrixPhy <- matrix(read_table(here::here("inst", "extdata", "cntPA.txt"), col_names = FALSE )$`X1`, 25, 25)  # contact matrix for conversational contacts
# Both
contactMatrixBoth <- contactMatrixCon + contactMatrixPhy

# my age strat: c(0, 1/12, 2/12, 3/12, 4/12, 5/12, 6/12, 7/12, 8/12, 9/12, 10/12, 11/12, 1, 2, 3, 4, 5, 10, 15, 25, 35, 45, 55, 65, 75)
# Julien age strat: c(0-23 months, 2, 5, 10, 20, 40, 60, 65, 75)

# GENERATE NEW MATRIX
# Change A: for 1yr-> 13-24 months
aug_rows <- rbind(
    contactMatrixBoth[1:12, ],
    t(matrix(rep(contactMatrixBoth[13, ], 12), 25, 12)),
    contactMatrixBoth[14:25, ]
)
contactMatrixBothA <- cbind(
    aug_rows[, 1:12],
    t(matrix(rep(aug_rows[, 13] / 12, 12), 12, 36)),
    aug_rows[, 14:25]
)

# Change B: 2,3,4 -> 2–4 years
aug_rows <- rbind(
    contactMatrixBothA[1:24, ],
    contactMatrixBothA[25:27, ] %>% apply(2, mean) ,
    contactMatrixBothA[28:36, ]
)
contactMatrixBothB <- cbind(
    aug_rows[, 1:24],
    matrix(aug_rows[, 25:27] %>% apply(1, sum), 34, 1),
    aug_rows[, 28:36]
)

# 5-9 (they are the same!)
#contactMatrixBothB[26, ]
#contactMatrixBoth[17, ]

# Other age groups
aug_rows <- rbind(
    contactMatrixBothB[1:26, ],
    # 10 - 20
    (contactMatrixBothB[27, ] + contactMatrixBothB[28, ] / 2) / 1.5,
    # 20 - 40
    (contactMatrixBothB[28, ] / 2 + contactMatrixBothB[29, ] + contactMatrixBothB[30, ] / 2) / 2,
    # 40 - 60
    (contactMatrixBothB[30, ] / 2 + contactMatrixBothB[31, ] + contactMatrixBothB[32, ] / 2) / 2,
    # 60 - 64
    contactMatrixBothB[32, ],
    # 65 - 74, 75 + 
    contactMatrixBothB[33:34, ]
)

contactMatrixBothC <- cbind(
    aug_rows[, 1:26],
    # 10 - 20
    (aug_rows[, 27] + aug_rows[, 28] / 2),
    # 20 - 40
    (aug_rows[, 28] / 2 + aug_rows[, 29] + aug_rows[, 30] / 2),
    # 40 - 60
    (aug_rows[, 30] / 2 + aug_rows[, 31] + aug_rows[, 32] / 2),
    # 60 - 64
    aug_rows[, 32] / 2,
    # 65 - 74, 75 + 
    aug_rows[, 33:34]
)

# Save file
write.csv(contactMatrixBothC, here::here("inst", "extdata", "cntAllv2.csv"))