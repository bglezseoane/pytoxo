# -*- coding: utf-8 -*-

###########################################################
# PyToxo
#
# A Python library for calculating penetrance tables of any
# bivariate epistasis model.
#
# Copyright 2021 Borja González Seoane
#
# Contact: borja.gseoane@udc.es
###########################################################

"""Graphical user interface entry point."""

import PySimpleGUI as sg

# ####################### GUI DESIGN ######################
# General style settings
sg.theme("SystemDefaultForReal")
window_general_font = "Verdana, 13"

# Main menu
menu = sg.Menu(
    [
        ["File", ["Open", "Save", "Exit"]],
        ["Edit", ["Paste", "Undo"]],
        ["Help", "About PyToxo GUI"],
    ]
)

# Epistatic model table
headings = [
    " Genotype definition ",
    "Penetrance expression",
    "Calculated penetrance",
]
rows = [[col + row for col in range(len(headings))] for row in range(100)]

# Configuration frame
configuration_frame = sg.Frame(
    title="Configuration",
    layout=[
        [
            sg.Combo(
                ("Heritability", "Prevalence"),
                key="prev_or_her_cb",
                size=(10, 1),
                readonly=True,
            ),
            sg.InputText(key="prev_or_her_in", size=(5, 1)),
        ],
        [
            sg.Text("MAFs"),
            # sg.InputText(key="mafs_in", size=(5, 1)),
            sg.Text("None model loaded", text_color="grey"),
        ],
    ],
)

# Model frame
model_frame = sg.Frame(
    title="Epistatic model",
    layout=[
        # [
        #     sg.Table(
        #         headings=headings,
        #         values=rows,
        #         vertical_scroll_only=True,
        #         justification="center",
        #     )
        # ]
        [sg.Text("None model loaded", text_color="grey")]
    ],
)

# Layout composition
layout = [
    [menu],
    [model_frame],
    [configuration_frame],
    [sg.Button("Calculate table")],
]

# Window composition
window = sg.Window(
    "PyToxo GUI",
    layout,
    default_element_size=(40, 1),
    grab_anywhere=False,
    font=window_general_font,
)
# #########################################################

# GUI event loop
def main():
    while True:
        event, values = window.read()

        # if event == "SaveSettings":
        #     filename = sg.popup_get_file("Save Settings", save_as=True, no_window=True)
        #     window.SaveToDisk(filename)
        #     # save(values)
        # elif event == "LoadSettings":
        #     filename = sg.popup_get_file("Load Settings", no_window=True)
        #     window.LoadFromDisk(filename)
        #     # load(form)
        if event == "Resize table":
            filename = sg.popup_get_file("Save Settings", save_as=True, no_window=True)
            window.SaveToDisk(filename)
            # save(values)
        elif event in ("Exit", None):
            break

    window.close()


if __name__ == "__main__":
    main()
