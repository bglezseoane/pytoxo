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


def main():
    # Style
    sg.theme("SystemDefaultForReal")
    window_general_font = "Verdana, 13"

    # Main menu
    menu_def = [
        ["File", ["Open", "Save", "Exit"]],
        ["Edit", ["Paste", "Undo"]],
        ["Help", "About PyToxo GUI"],
    ]

    # Model table components
    headings = [
        " Genotype definition ",
        "Penetrance expression",
        "Calculated penetrance",
    ]
    rows = [[col+row for col in range(len(headings))] for row in range(100)]

    # Main layout
    layout = [
        [sg.Menu(menu_def)],
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
            sg.Text("MAFs (separated with commas)"),
            sg.InputText(key="mafs_in", size=(5, 1)),
        ],
        [
            sg.Text("Interaction order"),
            sg.InputText(key="order_in", size=(2, 1)),
            sg.Button("Resize table"),
        ],
        [
            sg.Button("Calculate table"),
        ],
        [sg.Table(headings=headings, values=rows, vertical_scroll_only=True,
                  justification="center")],
    ]

    window = sg.Window(
        "PyToxo GUI",
        layout,
        default_element_size=(40, 1),
        grab_anywhere=False,
        font=window_general_font,
    )

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
