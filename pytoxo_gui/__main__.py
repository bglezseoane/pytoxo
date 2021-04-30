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
    sg.theme("DarkBlue4")
    window_general_font = "Verdana, 20"

    # Model table components
    headings = [
        "Genotype",
        "Penetrance expression",
        "Calculated penetrance probability",
    ]
    header = [[sg.Text(h, size=(20, 2), pad=(0, 0)) for h in headings]]

    input_rows = [
        [sg.Input(size=(20, 1), pad=(0, 0)) for col in range(len(headings))]
        for row in range(10)
    ]

    # Main layout
    layout = [
        [
            sg.Combo(
                ("Heritability", "Prevalence"),
                key="prev_or_her_cb",
                size=(20, 1),
                readonly=True,
            ),
            sg.InputText(key="prev_or_her_in"),
            sg.Text("MAFs (separated with commas)"),
            sg.InputText(key="mafs_in"),
        ],
        [
            sg.Button("Calculate table"),
        ],
        [header + input_rows],
        [
            sg.CBox("Checkbox", key="cb1"),
            sg.CBox("My second checkbox!", key="cb2", default=True),
        ],
        [
            sg.Radio("My first Radio!     ", "RADIO1", key="rad1", default=True),
            sg.Radio("My second Radio!", "RADIO1", key="rad2"),
        ],
        [
            sg.MLine(
                default_text="This is the default Text should you decide not to type anything",
                size=(35, 3),
                key="multi1",
            ),
            sg.MLine(default_text="A second multi-line", size=(35, 3), key="multi2"),
        ],
        [
            sg.Combo(("Combobox 1", "Combobox 2"), key="combo", size=(20, 1)),
            sg.Slider(
                range=(1, 100),
                orientation="h",
                size=(34, 20),
                key="slide1",
                default_value=85,
            ),
        ],
        [
            sg.OptionMenu(
                ("Menu Option 1", "Menu Option 2", "Menu Option 3"), key="optionmenu"
            )
        ],
        [
            sg.Listbox(
                values=("Listbox 1", "Listbox 2", "Listbox 3"),
                size=(30, 3),
                key="listbox",
            ),
            sg.Slider(
                range=(1, 100),
                orientation="v",
                size=(5, 20),
                default_value=25,
                key="slide2",
            ),
            sg.Slider(
                range=(1, 100),
                orientation="v",
                size=(5, 20),
                default_value=75,
                key="slide3",
            ),
            sg.Slider(
                range=(1, 100),
                orientation="v",
                size=(5, 20),
                default_value=10,
                key="slide4",
            ),
        ],
        [sg.Text("_" * 80)],
        [sg.Text("Choose A Folder", size=(35, 1))],
        [
            sg.Text("Your Folder", size=(15, 1), justification="right"),
            sg.InputText("Default Folder", key="folder"),
            sg.FolderBrowse(),
        ],
        [
            sg.Button("Exit"),
            sg.Text(" " * 40),
            sg.Button("SaveSettings"),
            sg.Button("LoadSettings"),
        ],
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

        if event == "SaveSettings":
            filename = sg.popup_get_file("Save Settings", save_as=True, no_window=True)
            window.SaveToDisk(filename)
            # save(values)
        elif event == "LoadSettings":
            filename = sg.popup_get_file("Load Settings", no_window=True)
            window.LoadFromDisk(filename)
            # load(form)
        elif event in ("Exit", None):
            break

    window.close()


if __name__ == "__main__":
    main()
