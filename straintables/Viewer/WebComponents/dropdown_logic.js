 
function ToggleDropdown(region_index) {
    var dropdown = document.getElementById("region_" + region_index + "_dropdown");
    if (!dropdown.classList.contains('show')) { CloseAllDropdowns();}
    dropdown.classList.toggle("show");
}


function CloseAllDropdowns() {
    var cname = "dropdown-content";
    var dropdowns = document.getElementsByClassName(cname);
    var i;
    for (i = 0; i < dropdowns.length; i++) {
        var openDropdown = dropdowns[i];
        if (openDropdown.classList.contains('show')) {
            openDropdown.classList.remove('show');
        }
    }
}

window.onclick = function(event) {
    if (!event.target.matches('.regionbtn')) {
        console.log(".");
        CloseAllDropdowns();
    }
};
