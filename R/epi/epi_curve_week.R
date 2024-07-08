library(tidyverse)
library(janitor)
library(zoo)
library(ggthemes)

## Code used to create weekly epi curve (figure 2b)
dat_outbreak <- read_csv("data/clean/dat_animal.csv")
dat_human <- read_csv("data/clean/dat_human.csv")

## Only interested in confirmed and probable cases in Romblon
dat_outbreak <- dat_outbreak %>% filter(province == "Romblon") %>% filter(date <= as.Date("2023-09-30"))
dat_human <- dat_human %>% filter(province == "Romblon") %>% filter(date <= as.Date("2023-09-30"))


dat_human <- dat_human %>%
  mutate(date_bitten = ymd(when_bitten),
         date = date_died) ## plotting the date died rather than date bitten

### Combine animal and human case data within one data frame for epicurve
## Creating a "Human Death" case status so that the human case can be included on the same epi curve as the animal cases
dat_human <- dat_human %>%
  select(date, province,municipality,barangay,month_year) %>% 
  mutate(case_status = "Human Death") 

dat_outbreak <- dat_outbreak %>%
  select(barangay,municipality,province,case_status,date,month_year)

## Combining the data frames
dat_outbreak <- rbind(dat_outbreak,dat_human) %>%
  mutate(case_status_fct = 
           factor(case_status, 
                  levels = c("Human Death","Confirmed", "Probable")))

levels(dat_outbreak$case_status_fct) <- 
  c("Human Death", 'Confirmed (Animal)', "Probable (Animal)")

## Combining the data frames
dat_outbreak <- dat_outbreak %>% 
  mutate(case_type = ifelse(case_status == "Human Death", "Human", "Animal")) %>% 
  mutate(case_type = factor(case_type, 
                            levels = c("Human","Animal")))

levels(dat_outbreak$case_type) <- 
  c("Human ", "Animal")

plot_weekly_epi_curve <- function(dat_outbreak, dat_human){
  

  weekly_breaks <- seq.Date(
    from = floor_date(min(dat_outbreak$date-30, na.rm=T),   "week", week_start = 1), # monday before first case
    to   = floor_date(max(dat_outbreak$date) + 5, "week", week_start = 1), # monday after last case
    by   = "week")
  
  ggplot(data = dat_outbreak) +          # set data 
    geom_histogram(                      # add histogram
      mapping = aes(x = date, fill = case_status_fct),   # map date column to x-axis
      breaks = weekly_breaks,
      col = "white")+                   # cases binned every 7 days, starting from first case (!) 
    scale_x_date(
      expand            = c(0,0),           # remove excess x-axis space before and after case bars
      date_breaks       = "month",        # date labels and major vertical gridlines appear every 3 Monday weeks
      date_minor_breaks = "week",           # minor vertical lines appear every Monday week
      date_labels       = "%b %Y",
    ) +
    scale_fill_manual(values=c("#001219","#c1121f","#DE7002")) +
    labs(
      # title = "Weekly cases Romblon Province \n\n",ae2012
      y = "Number of cases",
      x = "",
      fill="") + 
    ylim(0,4) +
    geom_hline(yintercept = seq(1, 4, by = 1),colour = "white", linewidth = 0.3)  +
    theme_clean(base_size = 14) +
    theme(legend.position = "none", 
          legend.text = element_text(size=11),
          axis.title.y = element_text(size = 13),
          axis.text.x = element_text(size = 11),
          # axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1),
          plot.margin = margin(5,30,5,20), 
          panel.border = element_blank(),
          plot.background = element_blank())
}

plot_weekly_epi_curve(dat_outbreak,dat_human)

png(filename = "output/figures/fig2b.png",width = 800, height = 150)
plot_weekly_epi_curve(dat_outbreak,dat_human)
dev.off()


pdf(file = "output/figures/fig2b.pdf",width = 12, height = 2)
plot_weekly_epi_curve(dat_outbreak,dat_human)
dev.off()


