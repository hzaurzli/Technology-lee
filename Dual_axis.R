library(tidyverse)
library(ggplot2)
library(readxl)
library(tidyverse)

test <- read_excel("Desktop/2222.xlsx")

test1 <- gather(test, key = "type", value = "value", -group)

ggplot(data = test1, aes(x = group)) +
  geom_col(aes(y = value, fill = type), position = "dodge") +
  theme_bw()

## 
test2 <- mutate(test, Vin_trans = 60 * Vin ) %>%
  select(-Vin) %>%
  gather(key = "type", value = "value", -group) %>%
  mutate(group = factor(group, levels = c("leaf1", "leaf2", "leaf3", "leaf4", "root", "stem", "flower")))

ggplot(data = test2, aes(x = group)) +
  geom_col(aes(y = value, fill = type), position = "dodge") +
  scale_y_continuous(
    name = "expression",
    sec.axis = sec_axis(~./60,
                        name = "Vin")
  ) +
  theme_bw()
