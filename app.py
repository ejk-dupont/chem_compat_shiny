from shiny import App, reactive, render, ui
from shiny.types import ImgData
from shiny.types import FileInfo
import cirpy
import re
import pandas as pd
from rdkit.Chem import Draw
from rdkit import Chem
from pathlib import Path
import platform
import os
from utils import rg_list, rg_list_numsort, element_list, adsorbent_list, determinefunctionals, mixlist_add, mixlist_bulkadd, mixlist_remove, mixlist_moveup, mixlist_movedown, makechart, makeadsorbchart, rg_details, newrg_list, download_chart, otherdata_search

# A card component wrapper.
def ui_card(title, *args):
    return (
        ui.div(
            {"class": "card mb-4"},
            ui.div(title, class_="card-header"),
            ui.div({"class": "card-body"}, *args),
        ),
    )

app_ui = ui.page_fluid(
    ui.div(
        ui.h3("Reactive Chemicals: Chemical Compatibility Tool (v2.6)"),
        ui.p(""),
    ),
    ui.navset_tab(
        ui.nav_panel("🧪 Mixture Components",
            ui.div(
                ui.p(""),
                ui.p("""Use a CAS number (preferred) or material name to find the reactive functional groups.
                Or, enter your own custom chemical if one cannot be found. """),
                ui.p("When your mixture of components is complete, switch over to the Compatibility Chart tab to view the table."),
                ui.p(""),
                ui.h4("Add to your mixture components"),
            ),
                ui.navset_card_tab(
                    ui.nav_panel("🔎 Search NIH database",
                        ui.h5("Search NIH database"),
                        ui.panel_well(
                            {"style":"column-count:2"},
                            ui.column(12,
                                ui.input_text("ask", "Chemical Search", value = "67-64-1"),
                                {"style":"height:=300px"}
                            ),
                            ui.column(4,
                                ui.input_action_button("search","Search"),
                                {"style":"height:=300px","style":"padding-top:31px"}
                            ),
                        ),
                        ui.div(
                        ui.p(""),
                            ui.row(
                                ui.column(4,
                                    ui.input_select("name","Name",choices=[]),
                                    ui.p("CAS",{"style":"margin:4px"}), 
                                    ui.output_text_verbatim("cas", placeholder = True),
                                    ui.p("SMILES",{"style":"margin:4px"}), 
                                    ui.output_text_verbatim("smiles", placeholder=True),
                                    ui.input_selectize("functgroup","Reactive Groups", choices=rg_list(),multiple=True),
                                    ui.input_action_button("addtolist","Add to list", class_="btn-success"),
                                ),
                                ui.column(6, {"align":"center"},
                                    ui.output_image("structure"),
                                    ui.output_text("pka"),
                                ),
                            ),
                        ),
                    ),
                    ui.nav_panel("✏ Add custom chemical",
                        ui.h5("Add a custom chemical"),
                        ui.input_text("cust_name","Name*"),
                        ui.input_text("cust_cas","CAS"),
                        ui.input_text("cust_smiles","SMILES"),
                        ui.input_action_button("get_cust_smiles","Get Reactive Groups from SMILES"),
                        ui.input_selectize("cust_functgroup","Reactive Groups*", choices=rg_list(),multiple=True),
                        ui.input_action_button("cust_addtolist","Add to list", class_="btn-success"),
                    ),
                    ui.nav_panel("📝 Add list of chemicals",
                        ui.h5("Add a custom chemical"),
                        ui.p("Upload a chart excel file to quickly load all the chemicals in your mixture."),
                        ui.layout_sidebar( 
                            ui.panel_sidebar(
                                ui.input_file("file1", "Choose chart file (.xlsx)", accept=[".xlsx"], multiple=False),
                            ),
                            ui.panel_main(
                                ui.p("To properly upload a file,"),
                                ui.img({"src": "table_layout.png","height":"157px","width":"437px"}),
                                ui.tags.ul(
                                    ui.tags.li("File must include headers for Name, CAS, Reactive group numbers, and SMILES"),
                                    ui.tags.li("Only Name and Reactive Group Numbers are required"),
                                    ui.tags.li("The sheet must be named 'Mixture' "),
                                    ui.tags.li("If you have downloaded a chart from this tool previously, it is already setup to be uploaded."),
                                ),
                            ),
                        ),
                        ui.p(""),
                        ui.output_ui("debug_df"),
                        ui.p(""),
                        ui.input_action_button("bulk_addtolist","Add to list", class_="btn-success"),
                    ),
                ),
            ui.div(
                ui.p(""),
                ui.h4("Your mixture components"),
                ui.output_data_frame("mixture")
            ),
            ui.div(
                ui.p(""),
                ui.row(
                    ui.column(4,
                        ui.input_action_button("removefromlist","Remove from list", class_="btn-danger"),
                    ),
                    ui.column(1,
                        ui.input_action_button("moveuplist","⬆"),
                    ),
                    ui.column(1,
                        ui.input_action_button("movedownlist","⬇"),
                    ),
                ),
            ),
            ui.div(
                ui.p(""),
                ui.input_selectize("adsorbs","Add your adsorbents", choices=adsorbent_list(),multiple=True),
            ),
        ),
        ui.nav_panel("📅 Compatibility Chart", 
               ui.p(""),
               ui.h4("Compatibility Chart"),
               ui.p("The chemical compatibility chart provides conservative advice about the effects of mixing compounds together, and if they are compatible."),
               ui.download_button("downloadchart","Download"),
               ui.p(""),
               ui.output_table("chemchart"),
               ui.p(""),
               ui.p(""),
               ui.p(""),
               ui.output_table("adschart"),
               ui.p(""),
               ui.p(""),
               ui.p(""),
               ui.div(
                   ui.row(
                       ui.p("Get more information on a specific mix:"),
                       ui.column(4,
                            ui.input_select("stream1","Component 1",["placeholder"]),
                       ),
                       ui.column(4,
                            ui.input_select("stream2","Component 2",["placeholder"]),
                       ),
                   ),
                   ui.output_ui("stream_mix"),
                   ui.output_text_verbatim("comments"),
               ),
               #ui.output_data_frame("chemchart_comments")
               #ui.output_data_frame("chemchart"),
        ),
        ui.nav_panel("⚠ Reactive Groups", 
               ui.p(""),
               ui.h4("Reactive Groups"),
               ui.p("Find a description of the reactive groups. You can filter by reactive group categories as well."),
               ui.row(
                   ui.column(4,
                        ui.panel_well(
                            ui.input_select("catgselect","Select reactive group category",choices = element_list(),selected='All'),
                            ui.input_select("rgselect","",rg_list(),size=20),
                        ),
                   ),
                   ui.column(7,
                        ui.row(
                            ui.tags.h6(ui.output_text("rgselected")),
                            ui.p(""),
                               ui.output_image("rgimage"),
                                {"align":"center", "padding":"margin-bottom:0 !important"}
                        ),
                        ui.row(
                            ui_card("🔥 Flammability",
                                    ui.output_text("flaminfo"),
                            ),
                            ui_card("🛑 Reactivity",
                                    ui.output_text("reactinfo"),
                            ),
                            ui_card("⚠ Toxicity",
                                    ui.output_text("toxinfo"),
                            ),
                            ui_card("📝 Notes",
                                    ui.output_text("noteinfo"),
                            ),
                            ui_card("🧪 Examples",
                                    ui.output_text("exinfo"),
                            ),
                        ),
                    ),
               ),
        ),
        ui.nav_panel("❔ About and feedback",
            ui.p(""),
            ui.p("""Here are some tutorial slides if you want some extra guidance. I am also working on adding addition features. If you use this tool, please provide feedback on the 
                 current or desired features by filling out the short form below. All answers are anonymous."""),
            ui.tags.a("👩‍🏫 Tutorial Slides", href="https://dupont.sharepoint.com/:p:/r/teams/ReactiveChemicals/Shared%20Documents/Chemical%20Interaction%20Matrices%20Resources/TIP%20Chem%20Compat%20Shiny/24-01-08%202.4%20Chem%20Compat%20Chart%20How%20To.pptx?d=w1607444f178741a3893c0e94ceb6db85&csf=1&web=1&e=Z65zCT",target = "_blank"),
            ui.br("    "),
            ui.tags.a("📝 Feedback form", href="https://forms.office.com/r/ur83G0mLH3?origin=lprLink",target = "_blank"),
            ui.p(""),
            ui.p("""This tool is under development by the Reactive Chemicals team. 
                 Thank you for your patience while this program develops."""
            ),
            ui.p("Recent Updates to v2.6:"),
            ui.tags.ul(
                ui.tags.li("💻 pKa preditions: The tool uses an ", ui.tags.a("open pKa predictor based on graph-convolution neural networks",href="https://xundrug.cn/molgpka",target ="_blank"),
                           """. This uses a SMILES string to determine a molecular graph, in other words, by using the chemical structure, 
                           to predict pKa """, ui.tags.a("(https://pubs.acs.org/doi/10.1021/acs.jcim.1c00075)", href="https://pubs.acs.org/doi/10.1021/acs.jcim.1c00075",target ="_blank"), 
                           """. Even better, it allows for queries from programs, such as Python, making the integration into the program easy."""),
                ui.p(""),
                ui.tags.li("🔢 Oxidation states for metals: In the most recent version of RDkit, there is now a function get determine ",
                           ui.tags.a("oxidation number", href="https://rdkit.org/docs/source/rdkit.Chem.rdMolDescriptors.html#rdkit.Chem.rdMolDescriptors.CalcOxidationNumbers",target="_blank"),
                           """. To determine if a metal is “active” or not, the current oxidation state is compared to a list of “stable” oxidation states. 
                           If the oxidation state does not match , it is considered “active” and is then further labeled with the proper metal reactive group. 
                           If the oxidation state does match , it is considered stable and unactive. This way, fully oxidized and generally inert metals, such as aluminum oxide (Al2O3) no longer have the overly conservative label of active metals, 
                           while pure aluminum keep the reactive label."""),
                ui.p(""),
                ui.tags.li("""📇 External database of common materials: During a search, if the program receives no hits from the NIH database, it will search this external database. 
                           Common polymers, such as polypropylene, do not have SMILES strings as they have no defined structure, or no set of repeating units. However, these polymers are often found in DuPont. 
                           To remedy this, an external database of common chemicals not found in NIH but are often found in DuPont was created. This includes compounds such as agar, gasoline, polyester, tar … 
                           CAS numbers always work best for these searches."""),
            ),
            
            ui.p(" If you have any questions or comments, please contact Elsa Koninckx (elsa.j.koninckx@dupont.com)."),
        ),
    ),
)

def server(input, output, session):
    df = reactive.Value(pd.DataFrame(columns = ['Name','CAS','Reactive Groups','SMILES','Structure']))
    df_bulk = reactive.Value(pd.DataFrame(columns = ['Name','CAS','Reactive Group Numbers','SMILES']))
    df_chart = reactive.Value(pd.DataFrame())
    df_comments = reactive.Value(pd.DataFrame())
    smiles_str = reactive.Value("")
    casnum = reactive.Value("")
    rginfo = reactive.Value([None]*8)
    est_pka = reactive.Value("")

    @reactive.Effect #NIH search, updated name, CAS, SMILES, reactive groups
    @reactive.event(input.search, ignore_none=True)
    def _():
        names_list = cirpy.resolve(input.ask(),'names')
        if names_list == None:
            results_name, results_cas, rgs  = otherdata_search(input.ask())
            if rgs == "None":
                ui.notification_show(input.ask()+" was not found", duration=10,close_button=True)
                ui.update_select("name", choices = [])
                ui.update_selectize("functgroup",selected="")
            else:
                smiles_str.set("")
                names_list = [results_name]
                ui.update_select("name", choices = names_list, selected = results_name)
                casnum.set(results_cas)
                ui.update_selectize("functgroup",selected=rgs)
        else:
            if not re.match('[\d/-]+$', input.ask()):
                if isinstance(names_list,str):
                    names_list = [names_list]
                names_list.insert(0,input.ask())

            ui.update_select("name", choices = names_list, selected = names_list[0])
            
            smiles_str.set(cirpy.resolve(input.ask(),'smiles'))
            casnum.set(cirpy.resolve(input.ask(),'cas'))
            fgs, pka = determinefunctionals(smiles_str())
            est_pka.set(pka)
            ui.update_selectize("functgroup",selected=fgs)
    
    @output
    @render.text
    @reactive.event(input.search, ignore_none=True)
    def cas():
        if input.ask() in casnum():
            casnum.set(input.ask())
        elif casnum() is not None:
            if type(casnum()) is list:
                casnum.set(casnum()[0])
        else:
            casnum.set("")
        return casnum()
    
    @output
    @render.text
    @reactive.event(input.search, ignore_none=True)
    def smiles():
        if smiles_str() is not None:
            return smiles_str()
        
    @output
    @render.image
    @reactive.event(input.search)
    def structure():
        if (platform.architecture()[1] == 'WindowsPE'):
            file_path = "temp.png"
        else:
            file_path = str(Path(__file__).parent / "temp.png")

        Draw.MolToFile(Chem.MolFromSmiles(smiles_str()),file_path)
        img: ImgData = {"src": file_path,"height":"200px"}
        return img
    
    @output
    @render.text
    @reactive.event(input.search, ignore_none=True)
    def pka():
        if est_pka() == " ":
            return ("")
        else:
            return("estimated pKa(s): {}".format(est_pka()))
        
        
    @reactive.Effect
    @reactive.event(input.addtolist)
    def _():
        fgstr = "".join(input.functgroup())
        if fgstr == "":
            ui.notification_show("A functional group must be assigned", duration=10,close_button=True)
        elif casnum() == "":
            ui.notification_show("Search a chemical", duration=10,close_button=True)
        elif len(df().index) > 0 and df()['CAS'].eq(casnum()).any():
            ui.notification_show("This material is already in your mixture", duration=10,close_button=True)
        else:
            updated_df = mixlist_add(df(), input.name(), casnum(), input.functgroup(),smiles_str())
            df.set(updated_df)

    @reactive.Effect
    @reactive.event(input.get_cust_smiles)
    def _():
        try:
            fgs, pka = determinefunctionals(input.cust_smiles())
            ui.update_selectize("cust_functgroup",selected=fgs)
        except:
            ui.notification_show("No reactive groups found in SMILES", duration=10,close_button=True)

    @reactive.Effect
    @reactive.event(input.cust_addtolist)
    def _():
        cust_fgstr = "".join(input.cust_functgroup())
        if cust_fgstr == "":
            ui.notification_show("A functional group must be assigned", duration=10,close_button=True)
        elif input.cust_name() == "":
            ui.notification_show("A name must be assigned", duration=10,close_button=True)
        elif len(df().index) > 0 and df()['CAS'].eq(input.cust_cas()).any():
            ui.notification_show("This material is already in your mixture", duration=10,close_button=True)
        else:
            if input.cust_cas() == '':
               cust_cas = '-'
            else:
               cust_cas = input.cust_cas()
            if input.cust_smiles() == '':
                cust_smiles = '-'
            else:
                cust_smiles = input.cust_smiles()
            cust_updated_df = mixlist_add(df(), input.cust_name(), cust_cas, input.cust_functgroup(),cust_smiles)
            df.set(cust_updated_df)
    
    @output
    @render.ui
    def debug_df():
        if input.file1() is None:
            ui.notification_show("Please upload an excel file", duration=10,close_button=True)
        else:
            f: list[FileInfo] = input.file1()
            df_bulk.set(pd.read_excel(f[0]["datapath"],sheet_name="Mixture"))
            return ui.HTML(df_bulk().to_html(classes="table table-bordered"))
    
    @reactive.Effect
    @reactive.event(input.bulk_addtolist)
    def _():
            df_nulls = df_bulk().isnull().any()
            
            if df_nulls["Name"] == True:
                ui.notification_show("An entry is missing a name", duration=10,close_button=True)
            elif df_nulls["Reactive Group Numbers"] == True:
                ui.notification_show("An entry is missing Reactive Group Numbers", duration=10,close_button=True)
            else:
                df_updated_bulk = mixlist_bulkadd(df(),df_bulk())
                df.set(df_updated_bulk)

    @reactive.Effect
    @reactive.event(input.removefromlist)
    def _():
        try:
            index_i = input.mixture_selected_rows()[0]
            if input.mixture_selected_rows()[0] is not None:
                index_i = input.mixture_selected_rows()[0]
                updated_df2 = mixlist_remove(df(),index_i)
                df.set(updated_df2)
        except:
            ui.notification_show("Select a component to delete", duration=10,close_button=True)
    
    @reactive.Effect
    @reactive.event(input.moveuplist)
    def _():
        try:
            index_i = input.mixture_selected_rows()[0]
            if index_i > 0 and len(df().index) > 0:
                updated_df3 = mixlist_moveup(df(),index_i)
                df.set(updated_df3)
        except:
            ui.notification_show("Select a component to move", duration=10,close_button=True)

    @reactive.Effect
    @reactive.event(input.movedownlist)
    def _():
        try:
            index_i = input.mixture_selected_rows()[0]
            if index_i < len(df().index)-1 and len(df().index) > 1:
                updated_df4 = mixlist_movedown(df(),index_i)
                df.set(updated_df4)
        except:
            ui.notification_show("Select a component to move", duration=10,close_button=True)

    @output
    @render.data_frame
    def mixture() -> ui.Tag:
        return render.DataGrid(
            df(), row_selection_mode="single",
        )
    
#Compatibility Chart    
    @output
    #@render.data_frame
    @render.table
    def chemchart():
        chart_comments, chart, chart_style = makechart(df())
        df_comments.set(chart_comments)
        df_chart.set(chart)
        return chart_style
        #return render.DataGrid(chart)
    
    @output
    @render.table
    def adschart():
        if len(input.adsorbs()) > 0:
            ads_chart = makeadsorbchart(df(),input.adsorbs())
            return ads_chart
        else:
            return None
    
    @reactive.Effect
    def _():
        ui.update_select("stream1",label="Component 1",choices=df()['Name'],selected=None)

    @reactive.Effect
    def _():
        ui.update_select("stream2",label="Component 2",choices=df()['Name'],selected=None)

    @output
    @render.ui
    def stream_mix():
        if df_chart().size == 0:
            return ""
        else:
            if df_chart().iat[int(input.stream1()),int(input.stream2())] != "":
                output = df_chart().iat[int(input.stream1()),int(input.stream2())]
            else:
                output = df_chart().iat[int(input.stream2()),int(input.stream1())]
            if output == "X":
                return ui.p("X - not self reactive", style = "color:grey")
            elif output == "SR":
                return ui.p("SR - caution, may be self reactive", style = "color:#FFD700")
            elif output == "C":
                return ui.p("C - caution, may be reactive", style = "color:#FFD700")
            elif output == "Y":
                return ui.p("Y - compatible", style = "color:#228B22")
            elif output == "N":
                return ui.p("N - not compatible", style = "color:#FF6347")

    @output
    @render.text
    def comments():
        if df_comments().size == 0:
            return ""
        else:
            if df_comments().iat[int(input.stream1()),int(input.stream2())] != "":
                comment = df_comments().iat[int(input.stream1()),int(input.stream2())]
            else:
                comment = df_comments().iat[int(input.stream2()),int(input.stream1())]
            comment2 = [sub.replace(r'\n','\n') for sub in comment]
            return "".join(comment2)


# Reactive Groups
    @reactive.Effect
    @reactive.event(input.catgselect)
    def _():
        if input.catgselect() == "All - alphabetical":
           newrglist = rg_list()
        elif input.catgselect() == "All - numerical":
           newrglist = rg_list_numsort()
        else:
            newrglist = newrg_list(input.catgselect()) 
        ui.update_select("rgselect",label="",choices=newrglist,selected=newrglist[0])
    
    @reactive.Effect
    @reactive.event(input.rgselect)
    def _():
        rginfo.set(rg_details(input.rgselect()))
    
    @reactive.Effect
    @reactive.event(input.rgselect)
    def _():
        rginfo.set(rg_details(input.rgselect()))
    
    @output
    @render.text
    def rgselected():
        return input.rgselect()

    @output
    @render.image
    @reactive.event(input.rgselect)
    def rgimage():
        dir = Path(__file__).parent
        rg = str(input.rgselect())
        rgnum =  int(rg[rg.find("(")+1:rg.find(")")])
        rg_img: ImgData = {"src": str(dir / ("rgpics/{}.png").format(rgnum)),"height":"130px"}
        return rg_img

    @output
    @render.text
    def flaminfo():
        return rginfo()[2]
    
    @output
    @render.text
    def reactinfo():
        return rginfo()[3]
    
    @output
    @render.text
    def toxinfo():
        return rginfo()[4]
    
    @output
    @render.text
    def noteinfo():
        return rginfo()[5]

    @output
    @render.text
    def exinfo():
        return rginfo()[6]
    
    @render.download()
    def downloadchart():
        download_chart(df())
        return str( dir / "chart.xlsx")

dir = Path(__file__).parent
app = App(app_ui, server,static_assets=dir)

#/srv/shiny_apps/elsa_chem_compat