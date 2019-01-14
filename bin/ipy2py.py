#!/bin/bash

jupyter nbconvert --to python --stdout --template templates/script.tpl $1 | awk 'NF > 0 {blank=0} NF == 0 {blank++} blank < 2' > $2
